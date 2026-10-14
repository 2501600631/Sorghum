# Current author-confirmed methods: 14 October 2026

The current methods describe the 98 samples already measured. They correct the previous 10-mesh description to 40 mesh, establish the origin of each stored spectrum as a two-repacking mean. The four reference constituents are wet-basis mass percentages, each averaged over two independently weighed test portions. Reference values and saved model results are unchanged.

## Moisture

Approximately 5.0 g of mixed powder was spread evenly in the sample pan of a Sartorius MA35M-000230V1 infrared moisture analyser. The temperature was 105 degrees Celsius and the `End of analysis` value was 0.0 min, selecting fully automatic termination based on the drying mass-loss curve. Moisture percentage is 100 times the heating mass loss divided by the initial sample mass.

## Crude protein

The Kjeldahl procedure referred to GB 5009.5-2025. A 1.000-g ground, mixed portion was weighed to 0.001 g. Digestion used 0.4 g CuSO4*5H2O, 6.0 g K2SO4 and 20 mL concentrated H2SO4 (at least 98%). Heating was increased gradually to 420 degrees Celsius; digestion continued for at least 1 h until the solution was green, clear and transparent. After cooling, the K9840 performed alkaline distillation with 400 g/L NaOH. Ammonia was collected in 20 g/L boric acid and titrated with 0.1000 mol/L HCl. Mixed indicator B used methyl red and bromocresol green at a 1:5 volume ratio; the endpoint was a stable pale grey-red colour after green. A full-process reagent blank was included. Blank-corrected nitrogen percentage was multiplied by 6.25.

For sample mass m in grams, HCl concentration c in mol/L and sample/blank titration volumes V and Vb in mL:

```text
crude_protein_percent = (1.4 * (V - Vb) * c / m) * 6.25
```

## Crude ash

Approximately 5 g was placed in a crucible previously incinerated to constant mass. The sample was carbonised and then incinerated in the SX2-4-10 muffle furnace at 550 degrees Celsius for 4 h. After cooling in a desiccator and weighing, incineration was repeated until successive masses differed by at most 0.5 mg. Ash percentage is 100 times the ash residue mass divided by the initial test-portion mass.

## Crude fibre

The F800 acid-alkali procedure referred to ISO 6865:2000. Approximately 1.000 g passing a 1.0-mm aperture was digested with 150 mL of 0.130 mol/L H2SO4 and then 150 mL of 0.230 mol/L KOH, each gently boiling for 30 min. These are amount-of-substance concentrations, rather than mass fractions. Filtration used heat-resistant sintered crucibles with 40-100 micrometre pores and diatomaceous earth previously incinerated at 500 +/-25 degrees Celsius for 4 h. The residue was washed thoroughly, dried at 130 +/-2 degrees Celsius for at least 2 h, cooled in a desiccator and weighed. Drying was repeated to constant mass. Incineration at 500 +/-25 degrees Celsius continued for at least 2 h and was repeated to constant mass before cooling and weighing.

At least two full-process blanks per batch used the same crucible type, filter aid and processing steps. The sample dry-minus-ash mass difference was corrected by the mean corresponding blank difference:

```text
crude_fibre_percent = 100 * ((sample_dry_mass - sample_ash_mass)
                         - mean(blank_dry_mass - blank_ash_mass))
                         / initial_sample_mass
```

The dry and ash masses include the same crucible and filter aid. Drying the analytical residue does not change the denominator to a dry-basis test portion.

## Sample preparation and spectra

Sealed samples stored at -15 degrees Celsius were brought to 20-25 degrees Celsius until the sample reached room temperature and the outside package had no condensation. Mixed whole grains were measured first in the instrument-compatible rotating diffuse-reflectance cup. A representative portion was ground using a laboratory high-speed grinder. All retained material was reground until the entire portion passed a 40-mesh sieve with 0.425-mm openings, with no discarded fraction. Temperature rise and moisture loss were controlled. Powder was mixed, promptly sealed and measured after its temperature stabilised.

The sample covered the entire window at sufficient depth to prevent transmission through the layer, with similar filling height and packing across measurements. The surface was naturally levelled without deliberate compression, and cups were cleaned between samples. Instrument preparation included warm-up, wavelength/photometric checks and internal high-reflectance reference correction, with status checks between batches.

The FOSS NIRS DS2500 acquired 400-2499.5 nm diffuse-reflectance spectra at 0.5-nm sampling intervals (4200 variables). Each measurement averaged 32 internal subscans. Two independent repackings were acquired for each sample and physical form, then averaged wavelength by wavelength. Exported values were apparent absorbance A(lambda)=log10(1/R(lambda)), with dimensionless corrected relative diffuse reflectance R. The archived records had no SNV, MSC, smoothing or derivative treatment before export. The fixed model grid retains every fourth wavelength.

## Records and audit interpretation

Current structured records are reference_method_metadata_publication_20261014.csv/.json and spectral_acquisition_metadata_publication_20261014.json in experimental_data/. The original metadata templates and frozen audit receipts retain their historical hashes. Their empty protocol fields describe information available to the original computational run. Historical replicate-status receipts correctly describe the one available exported record per form; the current acquisition declaration explains how that record was produced.

The author confirmation supplies the method descriptions. Individual repacking spectra, original batch identifiers, duplicate assay values, titration/gravimetric/blank weights and quantitative batch QC measurements were outside the supplied modelling data. No such records were reconstructed. Individual provider assignments are left in their historical state, while aggregate provider information and the independence declaration are retained.
