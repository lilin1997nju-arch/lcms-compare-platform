# Offline Unimod snapshot

`unimod.xml` is an unmodified copy of the Unimod database distributed by
[Unimod](https://www.unimod.org/downloads.html), not a database authored by this
application. It is distributed with its original Design Science License in
`unimod-license.txt`. Copyright and warranty terms are those of that license
and the upstream source. The retrieval timestamp, original source URL and
SHA-256 are recorded in `unimod-snapshot.json`.

The application reads this local snapshot only; analysis does not contact
Unimod. A maintainer can explicitly refresh the original source and license
with `tools/update_unimod.py`. Filtering at search time does not modify the
original XML. Database approval is a property of a database record, not proof
that a modification exists in any analyzed sample.
