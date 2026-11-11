# kinf research drivers

Moved here on 2026-11-11 from the karyotype-inference repo (`~/Projects/karyotype-inference`, branch `fastkar-port`,
commit `e22426e`), which is now frozen. The Python scripts import the shipped modules in `inst/python/kinf/`, and the
R exports call the package's `prep_data_for_py()`. Run everything from the fastKar repo root with the kinf env:
`env -u PYTHONPATH <kinf python> bench/kinf/research/<script>.py ...` (the interpreter is given by `fastKar:::kinf_python()`).

| File | What |
|---|---|
| `export_v2.R` | SK-OV-3 4-chromosome region benchmark (33 truths x 4 depths + real data) |
| `export_v3.R` | 5 CCLE complex events (stages `cand`, `probe`, `export`) |
| `export_wg.R` | whole-genome SK-OV-3 at `pix.size` |
| `wg_search.py` | greedy / PT search on one problem folder with a truth (gap to the truth NLL, exact recovery) |
| `wg_analyze.py` | plots and summaries of `wg_search.py` results |
| `kspace.py` | estimates the number of distinct karyotypes of a problem from random wirings |
| `validate_delta.py` | incremental (delta) NLL vs full recompute over random moves |
| `validate_walkdist.py` | walk-traversal expected map vs Floyd-Warshall |
| `test_roundtrip.py`, `test_hic_model.py` | gwalks <-> permutation round trip; forward-model spot check |

Data: the problem folders and results stay in the karyotype-inference main checkout,
`~/Projects/karyotype-inference/temp_data/{v2,v3,wg_skov3_1000000}` (gitignored there; the export scripts write to that
root).

**These folders were exported before 2026-11-11 (tile objective) and have no `bins.csv`, so the current kinf refuses
to load them.** Re-export before use. The export scripts write into the existing folders, so first point `outdir` to a
new folder rather than overwriting the old exports.

Not moved (still in the frozen KI repo): `scripts/export_model.R`, which duplicates `scripts/export_kinf_model.R`;
`scripts/ref_nll.R` and `tests/validate_forward.py`, which use the old tile objective and are replaced by the
`test-kinf-python.R` test "kinf scores like fastKar".
