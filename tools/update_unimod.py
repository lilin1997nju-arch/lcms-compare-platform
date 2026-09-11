"""Maintainer-only download of the original offline Unimod XML and license.

The application never calls this script or accesses Unimod over the network.
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen
import xml.etree.ElementTree as ET


def main():
    destination = Path(__file__).resolve().parents[1] / "lcms_feature_mvp/core/data"
    source = "https://www.unimod.org/xml/unimod.xml"
    with urlopen(source, timeout=60) as response:
        data = response.read()
    root = ET.fromstring(data)
    count = len(root.findall("{*}modifications/{*}mod"))
    if count < 1000:
        raise ValueError("Incomplete Unimod download; existing snapshot was not changed")
    with urlopen("https://www.unimod.org/dsl.txt", timeout=60) as response:
        license_data = response.read()
    if b"DESIGN SCIENCE LICENSE" not in license_data.upper():
        raise ValueError("Unimod license download is invalid")
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "unimod.xml").write_bytes(data)
    (destination / "unimod-license.txt").write_bytes(license_data)
    metadata = dict(source=source, retrieved_at=datetime.now(timezone.utc).isoformat(),
                    sha256=hashlib.sha256(data).hexdigest(), record_count=count,
                    license="Design Science License", source_modified=False)
    (destination / "unimod-snapshot.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
