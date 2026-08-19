# Vedolizumab sequence provenance and search assumption

- Sequence source: WHO Drug Information, Recommended INN List 62, Vol. 23 No. 3 (2009), vedolizumab entry.
- Source URL: https://cdn.who.int/media/docs/default-source/international-nonproprietary-names-%28inn%29/rl62.pdf
- Mature heavy chain: 451 amino acids.
- Mature light chain: 219 amino acids.
- Cross-check: Japanese Accepted Names database, Vedolizumab (Genetical Recombination), CAS 943609-66-3.
- Cross-check URL: https://jpdb.nihs.go.jp/jan/DetailList_ja?keyword=Vedolizumab+%28Genetical+Recombination%29&submit=all_alp%E6%A4%9C%E7%B4%A2

## Analysis scope

The public Vedolizumab sequence is used as the theoretical target for both raw files. QL2519 has not been independently sequence-confirmed from a sponsor-controlled FASTA, so matches in that sample are similarity evidence against the Vedolizumab target, not formal confirmation that its complete primary sequence is identical.

## Cys alkylation model check

Two otherwise identical searches were compared:

- Fixed Carbamidomethyl@C (+57.021464 Da): 5,302 accepted PSMs; 1,988 Cys-containing PSMs; 3,314 non-Cys PSMs.
- No fixed Cys modification: 3,319 accepted PSMs; 16 Cys-containing PSMs; 3,303 non-Cys PSMs.

The near-identical non-Cys yield and the large selective gain for Cys-containing spectra support fixed Carbamidomethyl@C as the data-supported working model. This is an analytical inference because the sample-preparation record is unavailable.
