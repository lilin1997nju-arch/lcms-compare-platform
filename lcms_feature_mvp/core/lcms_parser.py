#!/usr/bin/env python3
"""LC-MS import helpers.

Thermo RAW files are proprietary and need a converter such as ProteoWizard
or ThermoRawFileParser before scan arrays can be read. This module therefore
supports simple centroid CSV directly and provides a deterministic mock scan
generator for downloaded RAW files until the converter layer is available.
"""

from __future__ import annotations

import csv
import base64
import hashlib
import json
import math
import pickle
import random
import shutil
import struct
import tempfile
import zlib
from pathlib import Path
from typing import Callable
from xml.etree import ElementTree as ET

from .lcms_models import LCMSRawFile, LCMSSpectrumScan


MZML_CACHE_VERSION = 2


def infer_sample_id(path: Path) -> str:
    stem = path.stem
    if stem.endswith("-metadata"):
        stem = stem[:-9]
    if "MabThera" in stem:
        product = "MabThera"
    elif "Reditux" in stem:
        product = "Reditux"
    else:
        return stem
    replicate = stem.rsplit("_", 1)[-1] if "_" in stem else "1"
    return f"{product}_{replicate}"


def cv_value(element: ET.Element, accession: str) -> str | None:
    for child in element.iter():
        if strip_ns(child.tag) == "cvParam" and child.attrib.get("accession") == accession:
            return child.attrib.get("value", "")
    return None


def strip_ns(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def decode_binary_array(binary_array: ET.Element) -> tuple[str | None, list[float]]:
    array_type: str | None = None
    compressed = False
    precision = "64"
    binary_text = ""
    for child in binary_array:
        name = strip_ns(child.tag)
        if name == "cvParam":
            accession = child.attrib.get("accession")
            if accession == "MS:1000514":
                array_type = "mz"
            elif accession == "MS:1000515":
                array_type = "intensity"
            elif accession == "MS:1000574":
                compressed = True
            elif accession == "MS:1000521":
                precision = "32"
            elif accession == "MS:1000523":
                precision = "64"
        elif name == "binary":
            binary_text = child.text or ""
    if not binary_text:
        return array_type, []
    payload = base64.b64decode(binary_text)
    if compressed:
        payload = zlib.decompress(payload)
    fmt = "<f" if precision == "32" else "<d"
    width = 4 if precision == "32" else 8
    if len(payload) % width:
        raise ValueError("mzML binary payload size is not aligned with declared precision")
    return array_type, [value[0] for value in struct.iter_unpack(fmt, payload)]


def _mzml_cache_dir(path: Path) -> Path:
    """Return the generated, side-by-side cache directory for one mzML file."""
    return path.with_name(f".{path.name}.lcms-cache")


def _mzml_cache_signature(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


def _cached_mzml_levels(
    path: Path,
    project_id: str,
    ms_levels: tuple[int, ...] | None,
) -> tuple[LCMSRawFile, list[LCMSSpectrumScan]] | None:
    """Load selected scans from a valid side-by-side cache, if available.

    The cache stores one pickle stream per MS level.  Keeping the levels in
    separate streams lets the MS1 and MS2 stages read only the arrays they
    need, while the first stage still creates both streams in one XML pass.
    """
    cache_dir = _mzml_cache_dir(path)
    manifest_path = cache_dir / "manifest.json"
    if not manifest_path.exists():
        return None
    # The level-specific streams preserve XML order within one level. Keep
    # multi-level reads on the parser path so the public read_mzml ordering
    # remains unchanged for callers that request mixed levels.
    if ms_levels is None or len(set(ms_levels)) != 1:
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest.get("version") != MZML_CACHE_VERSION
            or manifest.get("source") != _mzml_cache_signature(path)
        ):
            return None
        level_files = dict(manifest.get("levels") or {})
        level_counts = {str(key): int(value) for key, value in dict(manifest.get("counts") or {}).items()}
        wanted = (
            sorted(int(level) for level in level_files)
            if ms_levels is None
            else sorted(set(int(level) for level in ms_levels))
        )
        if any(str(level) not in level_files for level in wanted):
            return None
        scans: list[LCMSSpectrumScan] = []
        for level in wanted:
            cache_path = cache_dir / str(level_files[str(level)])
            if not cache_path.is_file():
                return None
            loaded_count = 0
            with cache_path.open("rb") as handle:
                while True:
                    try:
                        scan = pickle.load(handle)
                    except EOFError:
                        break
                    if not isinstance(scan, LCMSSpectrumScan):
                        return None
                    scans.append(scan)
                    loaded_count += 1
            if str(level) in level_counts and loaded_count != level_counts[str(level)]:
                return None
        sample_id = infer_sample_id(path)
        return summarize_raw_file(path, sample_id, project_id, scans, parser_status="mzML-cache"), scans
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError, pickle.PickleError, EOFError):
        # A partial or stale cache is never fatal; the caller reparses mzML.
        return None


def _new_mzml_cache_dir(path: Path) -> Path:
    cache_dir = _mzml_cache_dir(path)
    return Path(tempfile.mkdtemp(prefix=f"{cache_dir.name}.", dir=str(path.parent)))


def read_mzml(
    path: Path,
    project_id: str = "lcms_mvp",
    ms_levels: tuple[int, ...] | None = (1,),
    use_cache: bool = True,
) -> tuple[LCMSRawFile, list[LCMSSpectrumScan]]:
    """Read selected centroid MS levels from mzML generated by ThermoRawFileParser."""
    path = Path(path)
    if use_cache:
        cached = _cached_mzml_levels(path, project_id, ms_levels)
        if cached is not None:
            return cached
    sample_id = infer_sample_id(path)
    raw_file_id = path.stem
    scans: list[LCMSSpectrumScan] = []
    selected_levels = set(ms_levels) if ms_levels is not None else None
    cache_dir: Path | None = _new_mzml_cache_dir(path) if use_cache else None
    cache_handles: dict[int, object] = {}
    cached_levels: set[int] = set()
    cached_counts: dict[int, int] = {}
    context = ET.iterparse(path, events=("end",))
    try:
        for _, element in context:
            if strip_ns(element.tag) != "spectrum":
                continue
            ms_level_text = cv_value(element, "MS:1000511")
            ms_level = int(float(ms_level_text)) if ms_level_text not in (None, "") else 1
            rt_text = cv_value(element, "MS:1000016")
            if rt_text in (None, ""):
                element.clear()
                continue
            rt = float(rt_text)
            tic_text = cv_value(element, "MS:1000285")
            base_mz_text = cv_value(element, "MS:1000504")
            base_int_text = cv_value(element, "MS:1000505")
            arrays: dict[str, list[float]] = {}
            for binary_array in element.iter():
                if strip_ns(binary_array.tag) != "binaryDataArray":
                    continue
                array_type, values = decode_binary_array(binary_array)
                if array_type:
                    arrays[array_type] = values
            mz_array = arrays.get("mz", [])
            intensity_array = arrays.get("intensity", [])
            if len(mz_array) != len(intensity_array):
                raise ValueError(f"{path.name}: mz/intensity array length mismatch at RT {rt}")
            precursor = next((item for item in element.iter() if strip_ns(item.tag) == "precursor"), None)
            activation_method = next(
                (
                    method
                    for accession, method in (
                        ("MS:1000422", "HCD"),
                        ("MS:1000133", "CID"),
                        ("MS:1000598", "ETD"),
                        ("MS:1000250", "ECD"),
                    )
                    if cv_value(element, accession) is not None
                ),
                None,
            )
            precursor_mz_text = cv_value(element, "MS:1000744")
            precursor_charge_text = cv_value(element, "MS:1000041")
            precursor_intensity_text = cv_value(element, "MS:1000042")
            isolation_lower_text = cv_value(element, "MS:1000828")
            isolation_upper_text = cv_value(element, "MS:1000829")
            collision_energy_text = cv_value(element, "MS:1000045")
            scan = LCMSSpectrumScan(
                scan_id=element.attrib.get("id", f"{raw_file_id}_scan_{len(scans) + 1}"),
                raw_file_id=raw_file_id,
                sample_id=sample_id,
                rt=rt,
                ms_level=ms_level,
                mz_array=mz_array,
                intensity_array=intensity_array,
                tic=float(tic_text) if tic_text not in (None, "") else sum(intensity_array),
                base_peak_mz=float(base_mz_text) if base_mz_text not in (None, "") else None,
                base_peak_intensity=float(base_int_text) if base_int_text not in (None, "") else (max(intensity_array) if intensity_array else 0.0),
                precursor_scan_id=precursor.attrib.get("spectrumRef") if precursor is not None else None,
                precursor_mz=float(precursor_mz_text) if precursor_mz_text not in (None, "") else None,
                precursor_charge=int(float(precursor_charge_text)) if precursor_charge_text not in (None, "") else None,
                precursor_intensity=float(precursor_intensity_text) if precursor_intensity_text not in (None, "") else None,
                isolation_window_lower_offset=float(isolation_lower_text) if isolation_lower_text not in (None, "") else None,
                isolation_window_upper_offset=float(isolation_upper_text) if isolation_upper_text not in (None, "") else None,
                activation_method=activation_method,
                collision_energy=float(collision_energy_text) if collision_energy_text not in (None, "") else None,
            )
            if cache_dir is not None:
                if ms_level not in cache_handles:
                    cache_path = cache_dir / f"ms{ms_level}.pkl"
                    cache_handles[ms_level] = cache_path.open("wb")
                pickle.dump(scan, cache_handles[ms_level], protocol=pickle.HIGHEST_PROTOCOL)
                cached_levels.add(ms_level)
                cached_counts[ms_level] = cached_counts.get(ms_level, 0) + 1
            if selected_levels is None or ms_level in selected_levels:
                scans.append(scan)
            element.clear()
        if cache_dir is not None:
            for handle in cache_handles.values():
                handle.close()
            manifest = {
                "version": MZML_CACHE_VERSION,
                "source": _mzml_cache_signature(path),
                "levels": {str(level): f"ms{level}.pkl" for level in sorted(cached_levels)},
                "counts": {str(level): cached_counts[level] for level in sorted(cached_levels)},
            }
            (cache_dir / "manifest.json").write_text(
                json.dumps(manifest, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            target_dir = _mzml_cache_dir(path)
            if target_dir.exists():
                shutil.rmtree(target_dir, ignore_errors=True)
            cache_dir.replace(target_dir)
            cache_dir = None
    except Exception:
        for handle in cache_handles.values():
            try:
                handle.close()
            except OSError:
                pass
        raise
    finally:
        if cache_dir is not None and cache_dir.exists():
            shutil.rmtree(cache_dir, ignore_errors=True)
    return summarize_raw_file(path, sample_id, project_id, scans, parser_status="mzML"), scans


def summarize_raw_file(
    path: Path,
    sample_id: str,
    project_id: str,
    scans: list[LCMSSpectrumScan],
    parser_status: str = "parsed",
) -> LCMSRawFile:
    mz_values = [mz for scan in scans for mz in scan.mz_array]
    rt_values = [scan.rt for scan in scans]
    return LCMSRawFile(
        raw_file_id=path.stem,
        sample_id=sample_id,
        project_id=project_id,
        file_name=path.name,
        file_path=str(path.resolve()),
        data_format=path.suffix.lower().lstrip(".") or "unknown",
        mz_min=min(mz_values) if mz_values else 0.0,
        mz_max=max(mz_values) if mz_values else 0.0,
        rt_min=min(rt_values) if rt_values else 0.0,
        rt_max=max(rt_values) if rt_values else 0.0,
        scan_count=len(scans),
        parser_status=parser_status,
    )


def read_centroid_csv(path: Path, project_id: str = "lcms_mvp") -> tuple[LCMSRawFile, list[LCMSSpectrumScan]]:
    """Read long-form centroid CSV: scan_id, sample_id, rt, mz, intensity."""
    grouped: dict[tuple[str, str, float], list[tuple[float, float]]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"scan_id", "sample_id", "rt", "mz", "intensity"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path.name}: missing required columns: {sorted(missing)}")
        for row in reader:
            sample_id = row["sample_id"]
            scan_id = row["scan_id"]
            rt = float(row["rt"])
            grouped.setdefault((sample_id, scan_id, rt), []).append((float(row["mz"]), float(row["intensity"])))

    scans: list[LCMSSpectrumScan] = []
    sample_id = infer_sample_id(path)
    for (row_sample_id, scan_id, rt), pairs in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][2])):
        mz_array = [mz for mz, _ in pairs]
        intensity_array = [intensity for _, intensity in pairs]
        tic = sum(intensity_array)
        max_idx = max(range(len(intensity_array)), key=intensity_array.__getitem__) if intensity_array else -1
        scans.append(
            LCMSSpectrumScan(
                scan_id=scan_id,
                raw_file_id=path.stem,
                sample_id=row_sample_id,
                rt=rt,
                ms_level=1,
                mz_array=mz_array,
                intensity_array=intensity_array,
                tic=tic,
                base_peak_mz=mz_array[max_idx] if max_idx >= 0 else None,
                base_peak_intensity=intensity_array[max_idx] if max_idx >= 0 else 0.0,
            )
        )
        sample_id = row_sample_id
    return summarize_raw_file(path, sample_id, project_id, scans), scans


def gaussian(rt: float, center: float, width: float, height: float) -> float:
    return height * math.exp(-0.5 * ((rt - center) / max(width, 1e-9)) ** 2)


def deterministic_seed(path: Path) -> int:
    digest = hashlib.sha256(path.name.encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def mock_scans_from_raw(path: Path, project_id: str = "zenodo_5005513") -> tuple[LCMSRawFile, list[LCMSSpectrumScan]]:
    """Generate reproducible centroid MS1 scans from a vendor RAW placeholder.

    The downloaded file remains the traceable source artifact. Generated scans
    are intentionally marked as mock because no RAW conversion tool is present.
    """
    rng = random.Random(deterministic_seed(path))
    sample_id = infer_sample_id(path)
    is_reditux = sample_id.startswith("Reditux")
    replicate_jitter = (deterministic_seed(path) % 17 - 8) * 0.002
    features = [
        # mz, rt, width, MabThera height, Reditux height
        (548.3124, 4.80, 0.10, 2.8e5, 2.7e5),
        (732.4411, 5.35, 0.13, 1.8e5, 3.4e5),  # increased in Reditux
        (914.5567, 5.90, 0.11, 2.5e5, 7.0e4),  # decreased in Reditux
        (1102.6589, 6.55, 0.16, 1.2e5, 1.3e5),
        (1264.7722, 7.10, 0.12, 0.0, 1.4e5),  # new in Reditux
        (1430.8831, 7.65, 0.15, 9.0e4, 0.0),  # missing in Reditux
        (1625.9965, 8.25, 0.18, 7.0e4, 7.5e4),
    ]
    scans: list[LCMSSpectrumScan] = []
    for index in range(260):
        rt = 3.5 + index * 0.025
        mz_array: list[float] = []
        intensity_array: list[float] = []
        for mz, center, width, mab_height, red_height in features:
            height = red_height if is_reditux else mab_height
            if height <= 0:
                continue
            local_rt = center + replicate_jitter + (0.035 if is_reditux and abs(mz - 732.4411) < 0.01 else 0.0)
            intensity = gaussian(rt, local_rt, width, height) * rng.uniform(0.94, 1.06)
            if intensity > 300:
                mz_error = mz * rng.uniform(-4, 4) / 1_000_000.0
                mz_array.append(mz + mz_error)
                intensity_array.append(intensity)
        for _ in range(rng.randint(8, 16)):
            mz_array.append(rng.uniform(200, 2000))
            intensity_array.append(rng.uniform(20, 900))
        pairs = sorted(zip(mz_array, intensity_array), key=lambda item: item[0])
        mz_array = [mz for mz, _ in pairs]
        intensity_array = [intensity for _, intensity in pairs]
        tic = sum(intensity_array)
        max_idx = max(range(len(intensity_array)), key=intensity_array.__getitem__) if intensity_array else -1
        scans.append(
            LCMSSpectrumScan(
                scan_id=f"{path.stem}_scan_{index + 1:04d}",
                raw_file_id=path.stem,
                sample_id=sample_id,
                rt=rt,
                ms_level=1,
                mz_array=mz_array,
                intensity_array=intensity_array,
                tic=tic,
                base_peak_mz=mz_array[max_idx] if max_idx >= 0 else None,
                base_peak_intensity=intensity_array[max_idx] if max_idx >= 0 else 0.0,
            )
        )
    return summarize_raw_file(path, sample_id, project_id, scans, parser_status="mock_from_vendor_raw"), scans


def load_lcms_file(
    path: Path,
    project_id: str = "lcms_mvp",
    ms_levels: tuple[int, ...] | None = (1,),
) -> tuple[LCMSRawFile, list[LCMSSpectrumScan]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return read_centroid_csv(path, project_id)
    if suffix == ".mzml":
        return read_mzml(path, project_id, ms_levels=ms_levels)
    if suffix == ".raw":
        return mock_scans_from_raw(path, project_id)
    raise ValueError(f"Unsupported LC-MS input format for MVP: {path}")


def load_lcms_directory(
    input_dir: Path,
    project_id: str = "lcms_mvp",
    ms_levels: tuple[int, ...] | None = (1,),
    progress_callback: Callable[[int, int, Path], None] | None = None,
) -> tuple[list[LCMSRawFile], dict[str, list[LCMSSpectrumScan]]]:
    files = sorted(
        path for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".csv", ".raw", ".mzml"}
    )
    if not files:
        raise ValueError(f"No .raw or .csv LC-MS files found under {input_dir}")
    raw_files: list[LCMSRawFile] = []
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for index, path in enumerate(files, start=1):
        raw_file, scans = load_lcms_file(path, project_id, ms_levels=ms_levels)
        raw_files.append(raw_file)
        scans_by_sample[raw_file.sample_id] = scans
        if progress_callback is not None:
            progress_callback(index, len(files), path)
    return raw_files, scans_by_sample


def write_mock_centroid_csv(path: Path, raw_inputs: list[Path]) -> None:
    """Persist generated scan points for inspection and repeatable tests."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scan_id", "sample_id", "rt", "mz", "intensity"])
        for raw_path in raw_inputs:
            _, scans = mock_scans_from_raw(raw_path)
            for scan in scans:
                for mz, intensity in zip(scan.mz_array, scan.intensity_array):
                    writer.writerow([scan.scan_id, scan.sample_id, f"{scan.rt:.5f}", f"{mz:.6f}", f"{intensity:.6f}"])


def raw_file_payload(raw_file: LCMSRawFile) -> dict[str, object]:
    return json.loads(json.dumps(raw_file.__dict__, ensure_ascii=False))
