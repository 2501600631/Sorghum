# Sorghum near-infrared spectroscopy: 98 paired samples

This upload package accompanies **Effects of grain grinding on near infrared prediction of sorghum composition**. It supplies paired whole-grain and ground-sample spectra, four wet-basis reference constituents, the complete analysis code, frozen repeated nested-validation results and three figure source packages.

Dataset title: **Paired sorghum spectra, reference composition and frozen nested-validation results for 98 samples**.
Package version: **13.2.0-98-publication-2026-10-14-v3**.
Scientific framework: **13.2.0-98**.
Repository destination: https://github.com/2501600631/Sorghum. This version is prepared locally for upload.

## Study and measurements

The analytical units are 98 independently coded sorghum samples. Sichuan University of Science and Engineering supplied 20 brewing sorghum samples, Liaoning Academy of Agricultural Sciences supplied 20 and Shanxi Academy of Agricultural Sciences supplied 58. Full sample identifiers, including their final digit, are preserved. Samples were stored in sealed bags at -15 degrees Celsius.

The authors corrected the original experimental description on 14 October 2026. All ground portions passed a 40-mesh sieve with 0.425-mm openings; retained material was reground until the entire portion passed. Every archived spectrum is the wavelength-wise arithmetic mean of two independent repackings, each measured using 32 internal subscans. The legacy `replicate=1` field identifies that one exported mean spectrum. The individual repacking spectra and batch identifiers are outside the supplied modelling export.

Measurements used a FOSS NIRS DS2500 in diffuse reflectance. The 4200 exported wavelengths span 400-2499.5 nm at 0.5-nm intervals. Apparent absorbance is `log10(1/R)`, with dimensionless corrected relative reflectance `R`. The model loader selects every fourth wavelength, giving 1050 variables from 400 to 2498 nm at 2-nm intervals. The archive retains the original, unpreprocessed exported values.

Moisture, crude protein, crude ash and crude fibre are mass percentages on an as-is (wet) basis. Each reference concentration is the mean of two independently weighed portions, with 98 valid values per constituent. Moisture used a Sartorius MA35 at 105 degrees Celsius with fully automatic `End of analysis=0.0 min`. Crude protein used blank-corrected Kjeldahl nitrogen multiplied by 6.25, with procedures referring to GB 5009.5-2025. Crude ash used incineration gravimetry. Crude fibre used acid-alkali digestion referring to ISO 6865:2000 and full-process blank correction. The current method records give reagent concentrations, temperatures, timings and calculation denominators. The archive contains the final modelling reference values; individual titration, gravimetric, duplicate-assay and batch QC measurements were not supplied.

## Files

| Path | Contents |
| --- | --- |
| `raw_inputs/` | English workbooks with preserved sample IDs, numerical values and row order |
| `experimental_data/` | Prepared spectra and reference chemistry, alignment, historical templates and current author-confirmed assay/acquisition metadata |
| `experimental_code/` | Model fitting, nested selection, inference and selective-grinding analysis |
| `v13_2_results/` | The complete frozen results, including the 152400-row inner candidate trace |
| `figures/Figure_01/` to `figures/Figure_03/` | Source workbook, plotting program, caption, PNG, PDF and SVG for each figure |
| `provenance/` | Historical freeze inventories, publication transformations and current correction/preservation audit |
| `CITATION.cff` | Citation title, authors and package version |
| `DATA_DICTIONARY.md` | Variables, units and result interpretation |
| `PUBLICATION_METHODS_UPDATE_20261014.md` | Current assay and spectral-acquisition details |
| `VERIFYING_THE_RELEASE.md` | Verification commands and historical/current metadata distinction |
| `RIGHTS_AND_LICENSE_STATUS.md` | Current redistribution and licence status |
| `SHA256SUMS.txt` | SHA256 inventory of every package file except the inventory itself |

## Environment and verification

The formal run used Python 3.8.3, NumPy 1.23.5, pandas 1.5.3, SciPy 1.10.1, scikit-learn 1.1.3, matplotlib 3.7.2, openpyxl 3.1.5 and xlrd 2.0.2. Use Python 3.8-3.10 with the supplied requirements.

```text
python -m pip install -r requirements.txt
python verify_publication.py
```

The current publication verifier checks every file checksum, English-language content, code comments/docstrings, raw/prepared agreement, all 600 repeat metrics, nested splits, complete candidate coverage and the author-confirmed current methods. It prints a JSON report and leaves the release files unchanged. On Windows, `00_publication_check.bat` runs the same check.

Further checks can be run in a working copy:

```text
python verify_repository.py
python run_v13_1_validate.py --verify-checksums
python experimental_code/test_independent_samples.py
python run_v13_1_algorithm_selftest.py
```

The last two commands write check receipts in `docs/`. The `v13_1` filenames are compatibility entry points for framework 13.2.0-98. The historical strict submission-metadata command checks original, empty metadata templates and therefore remains unsuitable as the current publication declaration. Current records supplement the historical audit rather than altering it.

## Reproduce figures and analysis

Run plotting programs in a working copy because they replace their image exports:

```text
python figures/Figure_01/Figure_01.py
python figures/Figure_02/Figure_02.py
python figures/Figure_03/Figure_03.py
```

Each program validates its plotted calculations and exports PNG, PDF and SVG. The supplied images and numerical source workbooks preserve the verified manuscript figure package.

For a complete new model run:

```text
python reproduce_formal.py --output-dir rerun_results
```

The command creates a separate working copy, freezes its current code and inputs, executes the fixed formal analysis and exports results into a previously nonexistent directory. The full search is computationally expensive. The packaging revision checks the stored analysis and executable equivalence; it does not rerun that search. Numerical-library and platform differences can affect refit results, so compare metrics and sample-level predictions as well as hashes. The smoke entry point performs a reduced search. `prepare_data.py` regenerates prepared inputs and import receipts from the raw workbooks.

## Interpretation

The fixed design uses five outer folds repeated ten times, four inner folds and seed 20260830. The 58800 OOF predictions comprise four constituents, fifteen methods, 98 samples and ten repetitions. Samples are the observational units. Ground-assisted calibration uses both physical forms during training and whole-grain spectra for held-out prediction. Its usefulness differs across constituents; the frozen results include the deterioration for crude fibre.

Selective grinding is a separate retrospective routing analysis. Model disagreement ranks samples within each outer-test fold and constituent. At the requested 20% fraction, 20 of 98 samples use existing ground-sample predictions and 78 retain whole-grain predictions; all 98 are evaluated. There is no retraining on those 20 samples. Random routing has the same fold quotas. Descriptive curves average 500 random selections, whereas formal inference uses exact expected random squared error. The evidence is internal repeated nested validation of the development dataset; an independent external batch is outside this archive.

## Provenance and deposition

Raw and prepared numerical data, scientific model code and all frozen results are byte-identical to the preceding author-reviewed package. Current methods document the authors' corrections without changing targets, predictions or validation splits. Python comments and docstrings are absent; scientific definitions are documented in the manuscript and this repository. The existing 9.20 crude-protein value for sample 202512-06142 remains unchanged.

`provenance/Publication_Correction_Audit_20261014_v3.json` records the preserved files and changed publication material. Original freeze receipts retain their original hashes. `UPLOAD_MANIFEST_20261014.md` identifies the update relative to the verified repository snapshot and the full inner trace to include. The archive is ready for author-controlled upload; licences and repository publication are separate author actions.
