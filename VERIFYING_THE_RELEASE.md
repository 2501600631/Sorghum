# Verification and reproduction

## Verify the stored results

Run `python verify_release.py` from the repository root in the supplied Python environment. On Windows, use `00_verify_release.bat`. A successful run prints a JSON report and leaves repository files unchanged.

The check covers the exact SHA256 inventory, English-language content in text and workbook XML, Python syntax, absence of public code comments/docstrings, raw/prepared numerical equality, historical data/code bindings, all 600 saved repeat metrics, outer/inner split isolation, 152400 valid candidates, selective-grinding comparisons and current method/citation metadata. Git internals and Python bytecode caches are excluded from the package inventory.

The stored analysis contains 98 samples, ten outer repetitions, five outer folds, four inner folds, 58800 out-of-fold predictions, 200 selected-workflow records and 18000 representation-screening records. Recomputing metrics from stored predictions verifies the archive; fitting the full analysis is a separate operation.

## Run a new analysis

```text
python reproduce_formal.py --output-dir rerun_results
```

The destination must be new. This command creates a temporary working copy, freezes its current code and inputs, runs the full fixed search and exports new results. The 20 October 2026 documentation revision preserved the historical model fits. Platform and library differences can affect refits; compare predictions and performance summaries as well as file hashes.

Figure scripts regenerate and overwrite image exports. Run them in a working copy:

```text
python figures/Figure_01/Figure_01.py
python figures/Figure_02/Figure_02.py
python figures/Figure_03/Figure_03.py
```

## Optional checks and compatible entry points

| Command | Scope and files written |
| --- | --- |
| `python verify_publication.py` | Compatible alias for current release verification; prints a report |
| `python verify_repository.py` | Core integrity and stored-result checks; prints a report |
| `python run_v13_1_validate.py --verify-checksums` | Historical package checks; prints a report |
| `python experimental_code/test_independent_samples.py` | Independent-sample regression test; writes a receipt under `docs/` |
| `python run_v13_1_algorithm_selftest.py` | Algorithm-family implementation checks; writes a receipt under `docs/` |
| `python run_v13_1_smoke_test.py` | Reduced model search in `v13_2_smoke_results/` |
| `python prepare_data.py` | Regenerates prepared data and import receipts |

Run commands that write files in a working copy. An altered receipt or figure export changes its release checksum. The `v13_1` filenames are retained compatibility paths for framework 13.2.0-98.

## Historical provenance and current methods

Original run receipts and dated packaging records retain their original contents. `provenance/Publication_Transformation.json` maps original files to the English repository derivative. `provenance/Publication_Correction_Audit_20261014_v3.json` records the subsequent author-confirmed methods correction. `provenance/Repository_Update_20261020.json` records the current journal-independent documentation update.

The 14 October 2026 method records describe wet-basis composition, two assay portions, 40-mesh preparation and two repackings behind each exported mean spectrum. Their metadata version remains tied to that confirmation date. The current package version in `DATASET_METADATA.json` and `CITATION.cff` identifies this repository revision.

The historical strict submission-metadata command reads archived empty templates for source groups and protocol fields. Its missing-field result describes those historical inputs. Current method information is checked by `verify_release.py`. The historical formal-result `--formal` validator binds original code hashes; use it with the original frozen run or a new run carrying its own matching freeze receipt.

The scientific evidence consists of internal repeated nested validation and retrospective selective-grinding analysis. Individual laboratory measurements, repacking spectra and an independent external-validation batch are outside this archive. Current method records are author declarations, with their stated record availability preserved.
