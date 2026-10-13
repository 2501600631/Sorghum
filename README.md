# Sorghum near-infrared spectroscopy: 98 paired samples

This repository accompanies the LWT manuscript on whole-grain prediction, ground-sample prediction, ground-assisted calibration, and retrospective selective grinding. It contains the experimental inputs, reproducible analysis, frozen formal results, and the source data and Python program for each figure. It was prepared for the existing Sorghum repository; this local package has not been uploaded or assigned a DOI.

## Study and data

The analytical units are 98 independently coded sorghum samples. Each sample has one whole-grain spectrum, one spectrum from the ground sample, and four reference measurements: moisture, crude protein, crude ash, and crude fibre. Full sample identifiers are retained, including their final digit. Sichuan University of Science and Engineering supplied 20 brewing sorghum samples, Liaoning Academy of Agricultural Sciences supplied 20 samples, and Shanxi Academy of Agricultural Sciences supplied 58 samples; individual provider-to-sample assignments were not supplied and have not been inferred. Samples were stored in sealed bags at -15 degrees Celsius. Reference measurements used a rapid moisture analyser, the Kjeldahl method, a muffle furnace, and a fibre analyser, respectively. Contents are recorded in percent. Reporting basis and the nitrogen-to-protein conversion factor have not been supplied.

Input spectra span 400-2499.5 nm at 0.5-nm intervals (4200 wavelengths). The fixed modelling grid spans 400-2498 nm at 2-nm intervals (1050 wavelengths). The completed formal analysis used five outer folds repeated ten times, four inner folds, and the fixed seed 20260830. There are 58,800 saved out-of-fold predictions: four components, fifteen methods, 98 samples, and ten repeats.

## Contents

| Directory or file | Purpose |
| --- | --- |
| `raw_inputs/` | English-header workbooks with unchanged source sample IDs, numerical values, and row ordering |
| `experimental_data/` | Machine-readable paired spectra, reference chemistry, alignment, and original metadata fields |
| `experimental_code/` | Model fitting, nested model selection, statistical inference, utility analysis, and validation |
| `v13_2_results/` | Completed frozen formal results and their historical audit receipts |
| `figures/Figure_01/` to `figures/Figure_03/` | An Excel source-data workbook, a standalone Python program, detailed caption, and image/vector exports for each figure |
| `provenance/` | Original freeze hashes and the original-to-public transformation audit |
| `CITATION.cff` | Dataset/code citation metadata for the author-approved eight-author list; no publication DOI or public-release date is asserted |
| `DATA_DICTIONARY.md` | Variables, units, modelling branches, and result interpretation |
| `VERIFYING_THE_RELEASE.md` | What was checked and the distinction between the public derivative and the original frozen run |
| `RIGHTS_AND_LICENSE_STATUS.md` | Rights status; no unapproved licence is assigned |

## Environment

Use a separate Python 3.8-3.10 environment. The formal computation used Python 3.8.3, NumPy 1.23.5, pandas 1.5.3, SciPy 1.10.1, scikit-learn 1.1.3, matplotlib 3.7.2, openpyxl 3.1.5, and xlrd 2.0.2. The fixed requirements are included. A modern Python version may not support these historical dependency versions.

```text
python -m pip install -r requirements.txt
python verify_repository.py
python run_v13_1_validate.py --verify-checksums
python experimental_code/test_independent_samples.py
python run_v13_1_algorithm_selftest.py
```

The `v13_1` names are retained as English compatibility entry points; the actual framework version is 13.2.0-98. The first three commands check the public release. The independent-sample test writes its check receipt to `docs/`. Run tests in a working copy if you need to preserve an unchanged release checksum inventory.

## Reproduce figures

Each program locates the Excel workbook beside itself and validates the plotted calculations before exporting PNG, PDF, and SVG. It does not retrain models. The supplied PNG files are the exact embedded figures extracted from the author-revised manuscript. Figure reproduction can replace the supplied image files; use a working copy to retain the release checksums.

```text
python figures/Figure_01/Figure_01.py
python figures/Figure_02/Figure_02.py
python figures/Figure_03/Figure_03.py
```

## Reproduce the formal analysis

```text
python reproduce_formal.py --output-dir rerun_results
```

This creates a temporary working copy, freezes the public code and current inputs, performs the fixed formal analysis, and copies the newly generated results to a previously nonexistent output directory. It does not overwrite the supplied frozen results. The complete search is computationally expensive and was not repeated merely to translate filenames or remove comments. Minor numerical differences across operating systems or numerical-library builds are possible; compare metrics and sample-level predictions as well as file hashes.

The smoke-test entry point uses a reduced search and cannot establish the manuscript's numerical findings. `prepare_data.py` regenerates the prepared tables from the translated raw workbooks. Use it only in a working copy because it updates the import audit and prepared-data files.

## Interpretation

Ground-assisted calibration uses both forms in outer-training data; its outer-test predictions require whole-grain spectra only. Selective grinding is a separate retrospective simulation. Within each outer-test fold and for each component, the standard deviation across eight whole-grain model predictions supplies a disagreement score. Quotas are assigned within folds. The selected samples use existing predictions generated from their ground-sample spectra; unselected samples retain whole-grain predictions. At the requested 20% fraction, 20 samples are routed to ground predictions and 78 retain whole-grain predictions, and all 98 samples are evaluated. No model is retrained using those 20 samples. The random comparison has matched fold quotas. Its descriptive curve averages 500 random selections; inference uses the exact expected random squared error instead.

This is internal repeated nested cross-validation and retrospective analysis of a development dataset. It contains no independent external batch validation. The 98 samples, rather than the 58,800 predictions or the ten repeats, are the observational units. The public code has no comments or docstrings, as requested; model definitions and interpretation are documented here and in the manuscript.

## Provenance

All frozen numeric tables and modelling input CSV files are preserved byte-for-byte. Raw workbook headers and worksheet names are English; the reference workbook was exported from XLS to XLSX without changing data values. Comments and docstrings were removed from Python files. A stale test-only assertion for sample 202512-06142 was corrected to the actual frozen crude-protein value of 9.20. The model data and predictions were not corrected or deleted. `provenance/Publication_Transformation.json` records hashes and executable-AST checks. Historical run receipts remain historical and have not been re-signed with publication-derivative code hashes.

See `RIGHTS_AND_LICENSE_STATUS.md` before public redistribution. No publication, repository upload, licence grant, independent validation, or journal acceptance is implied by the existence of this package.
