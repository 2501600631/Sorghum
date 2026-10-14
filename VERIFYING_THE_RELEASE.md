# Verification of the current publication package

Run `python verify_publication.py` in the supplied Python environment. It checks the full SHA256 inventory, English-language content in plain text, gzip tables and all workbook XML, absence of Python comments/docstrings, raw/prepared numerical equality, frozen-result integrity, all 600 repeat metrics, outer/inner split isolation, 152400 valid candidate records and current assay/acquisition metadata. It prints its JSON receipt and leaves package files unchanged.

The current metadata record is an author declaration dated 14 October 2026. Verification checks that declaration and its consistency with the supplied modelling export; it does not reconstruct original instrument logs or individual laboratory measurements.

## Preserved computation

The formal study contains 98 independent samples, ten outer repetitions, five outer folds, 58800 OOF predictions, 200 selected-workflow records, 18000 representation-screening records and 152400 valid inner candidates. The original numerical data, model code, split records, result tables and source figure workbooks are preserved. Current methods specify wet-basis composition and the averaging of two assay portions, the 40-mesh powder preparation and two acquisition repackings behind each exported mean spectrum.

The packaging check recomputes saved performance metrics and inspects candidate/split records. The complete formal model search is available through reproduce_formal.py and has not been rerun during this documentation revision. The supplied design is internal repeated nested validation; an independent external-validation batch is outside the package.

## Historical hashes and validators

Historical receipts retain the hashes of the actual original run. `provenance/Publication_Transformation.json` maps those hashes to the earlier English publication derivative and records executable-AST preservation. `provenance/Publication_Correction_Audit_20261014_v3.json` maps the preceding publication package to this correction. Comments and docstrings were already absent from the public scientific code, and scientific executable behavior is unchanged in this revision.

Original data/code freeze inventories are retained in provenance/. Current derivative freeze inventories remain valid because the original computational files are unchanged. The new publication metadata and checker are covered by SHA256SUMS.txt, which inventories every delivered file except itself.

The original strict submission-metadata entry point reads archived empty templates rather than the current author declaration. It therefore retains its historical failure for missing per-sample source groups and original protocol fields. Use verify_publication.py for current package verification. The original strict formal-result validator also binds original code hashes, so its --formal mode is for an original frozen run or a new run carrying its own freeze receipt. Keeping that historical check preserves the original provenance; it is not changed into a claim that corrected documentation existed at the original run date.

## Optional executable checks

The independent-sample and twenty-family algorithm checks write receipts in docs/ and should be run in a working copy. Figure scripts validate their numerical source calculations and replace their PNG/PDF/SVG exports. Their prior successful check receipts are retained; each plotting source workbook and image remains part of the verified figure package.

Local verification prepares this archive for upload. Repository deposition and licence choice remain author-controlled actions recorded separately in PUBLICATION_METADATA.json and RIGHTS_AND_LICENSE_STATUS.md.
