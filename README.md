# Sorghum near-infrared spectroscopy: 98 paired samples

This repository provides paired whole-grain and ground-sample spectra, reference composition, analysis code, frozen repeated nested-validation results and three figure source packages for 98 sorghum samples. It supports reuse and verification independently of a submission venue.

- Dataset: **Paired sorghum spectra, reference composition and frozen nested-validation results for 98 samples**
- Repository package: **13.2.0-98-repository-2026-10-20**
- Scientific framework: **13.2.0-98**
- Repository: <https://github.com/2501600631/Sorghum>

## Study and measurements

The analytical units are 98 independently coded samples: 20 supplied by Sichuan University of Science and Engineering, 20 by Liaoning Academy of Agricultural Sciences and 58 by Shanxi Academy of Agricultural Sciences. Complete sample identifiers, including their final digit, are retained. Samples were stored in sealed bags at -15 degrees Celsius.

The FOSS NIRS DS2500 measured diffuse reflectance over 400-2499.5 nm at 0.5-nm intervals (4200 variables). Ground portions passed entirely through a 40-mesh sieve with 0.425-mm openings. Each archived spectrum is the wavelength-wise mean of two independent repackings, with 32 internal subscans per measurement. The legacy `replicate=1` field identifies the exported mean spectrum. Apparent absorbance is `log10(1/R)`, where `R` is dimensionless corrected relative reflectance. The model loader retains every fourth wavelength, producing 1050 variables at 2-nm intervals.

Moisture, crude protein, crude ash and crude fibre are mass percentages on an as-is (wet) basis. Each reference value averages two independently weighed portions; all four constituents have 98 valid values. The current author-confirmed assay and acquisition details are in [PUBLICATION_METHODS_UPDATE_20261014.md](PUBLICATION_METHODS_UPDATE_20261014.md). Those dated method records preserve their existing filenames and versions.

## Repository contents

| Path | Contents |
| --- | --- |
| `raw_inputs/` | Source workbooks preserving sample IDs, values and row order |
| `experimental_data/` | Prepared spectra and reference values, alignment, historical templates and dated method records |
| `experimental_code/` | Model fitting, nested selection, inference and selective-grinding analysis |
| `v13_2_results/` | Frozen results, including the complete 152400-row inner candidate trace |
| `figures/Figure_01/` to `figures/Figure_03/` | Source workbook, plotting program, caption, PNG, PDF and SVG |
| `provenance/` | Original freeze records and subsequent documentation audits |
| `DATASET_METADATA.json`, `CITATION.cff` | Current dataset description, authors and repository package version |
| `DATA_DICTIONARY.md` | Variables, units and result interpretation |
| `VERIFYING_THE_RELEASE.md` | Verification, reproduction and historical-entry-point guidance |
| `RIGHTS_AND_LICENSE_STATUS.md` | Licence and reuse status |
| `SHA256SUMS.txt` | Checksums for package files, excluding itself and Git internals |

## Verify the stored package

The historical formal run used Python 3.8.3, NumPy 1.23.5, pandas 1.5.3, SciPy 1.10.1, scikit-learn 1.1.3, matplotlib 3.7.2, openpyxl 3.1.5 and xlrd 2.0.2. Use Python 3.8-3.10 with the pinned requirements.

```text
python -m pip install -r requirements.txt
python verify_release.py
```

On Windows, `00_verify_release.bat` runs the same check. These are the stable verification entry points. The verifier checks file integrity, raw/prepared agreement, all 600 saved repeat metrics, outer/inner split isolation, candidate coverage, current method records and citation metadata. It prints a JSON report without fitting models or changing package files. The existing `verify_publication.py` and `00_publication_check.bat` remain compatible entry points.

## Reproduce figures or models

Run figure programs in a working copy because they replace their PNG, PDF and SVG exports:

```text
python figures/Figure_01/Figure_01.py
python figures/Figure_02/Figure_02.py
python figures/Figure_03/Figure_03.py
```

Run the complete fixed analysis into a new output directory:

```text
python reproduce_formal.py --output-dir rerun_results
```

This command creates a separate working copy and executes the full model search. It is computationally expensive. Compare sample-level predictions and metrics when numerical libraries or platforms differ. `prepare_data.py` regenerates prepared inputs and import receipts from the source workbooks. Additional checks and their output behaviour are described in [VERIFYING_THE_RELEASE.md](VERIFYING_THE_RELEASE.md).

## Evidence and interpretation

The fixed design uses five outer folds repeated ten times, four inner folds and seed 20260830. The 58800 out-of-fold predictions cover four constituents, fifteen methods, 98 samples and ten repetitions. Samples are the observational units. Ground-assisted calibration uses both forms during training and whole-grain spectra for held-out prediction. Its effect differs across constituents, including deterioration for crude fibre.

Selective grinding is a retrospective routing analysis. Model disagreement ranks samples within each outer-test fold and constituent. A requested 20% budget routes 20 of 98 samples to existing ground-sample predictions; the remaining 78 retain whole-grain predictions. All 98 are evaluated, with matched fold quotas for random routing. Descriptive curves average 500 random selections, while formal inference uses exact expected random squared error.

The evidence is internal repeated nested validation of the development dataset. The archive contains no independent external batch, individual repacking spectra, duplicate-assay measurements or quantitative batch QC records. Those evidence boundaries also apply to reuse of the frozen results.

## Version history and citation

The 20 October 2026 update generalises repository documentation, metadata and verification entry points. Raw inputs, prepared experimental records, scientific model code, figures, validation splits and frozen results retain their bytes from baseline commit `54d63d630d11daa65efcfc313b27192ee2b41947`; the full model search was not rerun.

Dates and the words `publication` or `submission` in older filenames, run profiles and receipts identify historical packaging or validation steps. Original records remain available under their existing paths. [provenance/Repository_Update_20261020.json](provenance/Repository_Update_20261020.json) records this documentation update; `PUBLICATION_METADATA.json` remains a compatibility copy of the current dataset metadata.

Use [CITATION.cff](CITATION.cff) to cite this dataset and analysis package. The archive has no assigned DOI or licence; see [RIGHTS_AND_LICENSE_STATUS.md](RIGHTS_AND_LICENSE_STATUS.md) for reuse terms.
