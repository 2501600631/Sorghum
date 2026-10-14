# Data dictionary

## Samples and reference chemistry

`raw_inputs/reference_chemistry_raw.xlsx` contains the `reference_chemistry` worksheet. Its columns are `sample_id`, `moisture_pct`, `crude_protein_pct`, `crude_ash_pct`, and `crude_fiber_pct`. The prepared table `experimental_data/reference_chemistry.csv` uses the corresponding English target keys without the `_pct` suffix. Values are reported in percent; prediction RMSE is in percentage points. The authors confirmed on 14 October 2026 that all four constituents are on an as-is (wet) basis; Kjeldahl crude protein uses N × 6.25. Each aggregate value is the mean of two separately weighed portions. See the new publication-method metadata for the confirmed details.

`sample_id` is the complete identifier of a unique analytical sample. A final 1 or 2 is part of that identifier. Neither suffix pairs nor matched whole/ground spectra are averaged into new samples. `sample_alignment.csv` gives the source Excel row for chemistry and each spectral form; Excel rows are one-based and include the header row.

## Spectra

`raw_inputs/paired_spectra_raw.xlsx` contains `whole_grain_400_2499p5` and `ground_grain_400_2499p5`. The first column is `sample_id`; numerical column headings give wavelengths in nm. There are 4200 spectral values per sample, from 400 to 2499.5 nm in steps of 0.5 nm. Spectral values are unpreprocessed apparent absorbance, A(lambda) = log10(1/R(lambda)), where R is dimensionless corrected diffuse reflectance relative to the instrument reference. Powder preparation used full passage through a 40-mesh sieve with 0.425-mm openings.

The two `experimental_data/*_grain_replicate_spectra.csv.gz` files contain `sample_id`, `replicate`, and columns named `wl_400.0` through `wl_2499.5`. The legacy field `replicate` is always 1 and identifies the single available exported mean spectrum in each form. The authors confirmed that this spectrum is the wavelength-wise arithmetic mean of two independent repackings, each obtained by averaging 32 internal subscans. Individual repacking spectra and acquisition-batch identifiers are outside the supplied modelling export. CSV gzip compression is lossless. The model loader selects every fourth exported wavelength, producing the fixed 2-nm modelling grid.

## Modelling branch keys

| Key | Meaning | Test input |
| --- | --- | --- |
| `nested_strong_whole` | Inner-selected model using whole-grain calibration spectra | Whole-grain spectrum |
| `ground_reference` | Independently inner-selected model using ground-sample calibration spectra | Ground-sample spectrum |
| `nested_ground_assisted` | Inner-selected EPO-PLS, EPO-elastic net, pooled elastic net, or diPLS using both forms during calibration | Whole-grain spectrum |
| `whole_pls`, `whole_ridge`, `whole_elastic`, `whole_ipls`, `whole_svr`, `whole_krr`, `whole_pls_krr`, `whole_extratrees` | Inner-selected predictions for the eight whole-grain model families | Whole-grain spectrum |
| `epo_pls`, `epo_elastic`, `pooled_elastic`, `dipls` | Inner-selected predictions for individual assisted families | Whole-grain spectrum |

## Out-of-fold prediction table

`v13_2_results/01_predictions/all_targets_outer_oof_predictions.csv` contains:

| Column | Meaning |
| --- | --- |
| `target` | One of `crude_protein`, `moisture`, `crude_ash`, and `crude_fiber` |
| `method` | Modelling branch or family key |
| `sample_id` | Complete held-out sample identifier |
| `inference_group` | Resampling group; the sample identifier when no verified source groups were available |
| `outer_repeat` | Zero-based outer repetition index, 0-9 |
| `outer_fold` | Zero-based outer test-fold index, 0-4 |
| `y_true` | Laboratory reference value, percent on an as-is (wet) basis |
| `y_pred` | Held-out model prediction, percent on an as-is (wet) basis |
| `selected_family` | Model family selected inside the outer-training data |
| `representation` | Wavelength-range and preprocessing key |
| `params_json` | Selected hyperparameters encoded as JSON |

Each component-method-sample combination has ten predictions. Within every repeat, each sample is held out once. `outer_repeat_performance.csv` computes metrics using the 98 held-out predictions in that repeat; manuscript performance tables report medians and interquartile ranges across ten such metrics. Averaging a sample's ten predictions is used for descriptive scatterplots and is not an increase in sample size.

## Selection and inference

`02_selection/` preserves training/test identifiers, inner splits, representation screening, selected workflows, and candidate tuning traces. `04_inference/` gives paired assisted-versus-whole comparisons. `05_utility_boundary/` gives descriptive selective-grinding curves, sample decisions, and oracle analyses. `06_robustness/` gives multiplicity-adjusted comparisons, practical-gain sensitivities, corrected repeated-CV sensitivity, and a fixed-prediction label-permutation diagnostic. A fixed-prediction permutation is not a permutation-based refit of the entire nested model search.

Selective-grinding budgets are requested fractions; realised fractions follow integer allocation. Requested 0.1 and 0.2 allocate 10/98 and 20/98 samples. Disagreement scores rank samples within each outer fold, separately by component and repeat. Positive `delta_mse_guided_better_positive` values indicate lower error for guided routing. Holm-adjusted P values control the family of eight component-by-budget comparisons. Confidence intervals and hypothesis tests follow the script definitions; descriptive random curves and exact expected random squared-error inference are different summaries.

## Figure workbooks

`Figure_01_data.xlsx` provides whole-grain OOF predictions, averaged sample predictions, and repeat metrics. `Figure_02_data.xlsx` provides workflow comparisons, paired inference, and plotted numerical summaries. `Figure_03_data.xlsx` provides routing curves, realised quotas, and selective-grinding comparisons. `ReadMe` or equivalent metadata sheets describe plotted variables and units. Numeric cell values are unchanged from the corresponding verified figure source-data package.

## Audit files and metadata

Blank `source_group`, cultivar, batch, reporting-basis, analytical-repeat and laboratory-protocol fields mean those fields were not provided in the frozen computational input. They have not been replaced with fabricated values. These historical template names are retained to preserve reproducibility. Confirmed aggregate sample-provider and reference-method information is stated in the README and manuscript. Original hashes identify the original frozen run; publication-derivative hashes identify this English release. Their different purposes are documented in `VERIFYING_THE_RELEASE.md`.

## Current and historical reference-method metadata

The original reference_method_metadata_template.csv and frozen audit records remain unchanged and reflect the inputs available during the historical formal run. They are not the current publication method declaration. Current author-confirmed information is in reference_method_metadata_publication_20261014.csv and its JSON companion. Current detailed digestion, titration and spectral-acquisition parameters are author-confirmed. The fibre result uses a full-process blank correction. Individual duplicate-assay values, raw titration or weighing records and batch QC numbers were not supplied. The available modelling values and frozen results are unchanged. Spectral acquisition is documented in spectral_acquisition_metadata_publication_20261014.json; historical technical_replicate_status.json describes the single exported modelling record, rather than the now-confirmed two acquisition repackings.
