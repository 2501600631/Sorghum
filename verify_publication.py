import csv
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import verify_repository as release

ROOT = Path(__file__).resolve().parent
VERSION = '13.2.0-98-publication-2026-10-14-v3'
TITLE = 'Paired sorghum spectra, reference composition and frozen nested-validation results for 98 samples'

def require(condition, message):
    if not condition:
        raise AssertionError(message)

def read_json(relative):
    return json.loads((ROOT / relative).read_text(encoding='utf-8'))

def check_inventory_coverage():
    listed = []
    for line in (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        if line:
            listed.append(line.split('  ', 1)[1])
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name != 'SHA256SUMS.txt'}
    require(len(listed) == len(set(listed)), 'Duplicate release inventory entries')
    require(set(listed) == actual, 'Release inventory does not cover the exact package contents')
    return len(actual)

def check_publication_methods():
    methods = read_json('experimental_data/reference_method_metadata_publication_20261014.json')
    require(len(methods) == 4, 'Current assay metadata row count')
    by_target = {row['target']: row for row in methods}
    require(set(by_target) == set(release.TARGETS), 'Current assay metadata target coverage')
    for row in methods:
        require(row['reporting_basis'] == 'as-is (wet)' and row['result_unit'] == 'mass percent', 'Reference basis or units')
        require(row['independent_test_portions'] == 2 and row['reference_aggregation'] == 'arithmetic mean', 'Reference aggregation')
        require(row['valid_reference_values'] == 98 and row['missing_reference_values'] == 0, 'Reference coverage')
        require(row['confirmation_date'] == '2026-10-14' and row['numerical_reference_values_changed'] is False, 'Confirmation provenance')
        require(row['individual_assay_records_in_archive'] is False, 'Raw laboratory record availability')
    moisture = by_target['moisture']
    require(moisture['instrument_model'] == 'MA35M-000230V1', 'Moisture analyser model')
    require(moisture['temperature_C'] == 105 and moisture['endpoint_value_min'] == 0.0 and moisture['endpoint_mode'] == 'fully automatic', 'Moisture operating settings')
    require(moisture['oven_method_validation_completed'] is False, 'Moisture metadata must not assert an unperformed comparison')
    protein = by_target['crude_protein']
    require(protein['reference_standard'] == 'GB 5009.5-2025' and protein['instrument_model'] == 'K9840', 'Protein method reference')
    for key, value in {'sample_mass_g': 1.000, 'weighing_precision_g': 0.001, 'copper_sulfate_pentahydrate_g': 0.4, 'potassium_sulfate_g': 6.0, 'concentrated_sulfuric_acid_mL': 20.0, 'digestion_temperature_C': 420, 'minimum_digestion_time_h': 1.0, 'sodium_hydroxide_g_L': 400, 'boric_acid_g_L': 20, 'titrant_mol_L': 0.1000, 'nitrogen_to_protein_factor': 6.25}.items():
        require(protein[key] == value, 'Protein metadata mismatch: ' + key)
    require(protein['full_process_reagent_blank'] is True and protein['methyl_red_to_bromocresol_green_volume_ratio'] == '1:5', 'Protein blank or indicator')
    ash = by_target['crude_ash']
    require(ash['incineration_temperature_C'] == 550 and ash['initial_incineration_time_h'] == 4.0 and ash['constant_mass_maximum_difference_mg'] == 0.5, 'Ash method settings')
    fiber = by_target['crude_fiber']
    require(fiber['reference_standard'] == 'ISO 6865:2000' and fiber['instrument_model'] == 'F800', 'Fibre method reference')
    require(fiber['concentration_basis'] == 'amount-of-substance concentration in mol/L', 'Fibre reagent concentration definition')
    for key, value in {'maximum_particle_sieve_aperture_mm': 1.0, 'sulfuric_acid_mol_L': 0.130, 'sulfuric_acid_mL': 150, 'potassium_hydroxide_mol_L': 0.230, 'potassium_hydroxide_mL': 150, 'acid_gentle_boiling_time_min': 30, 'alkali_gentle_boiling_time_min': 30, 'filter_aid_preincineration_temperature_C': 500, 'filter_aid_preincineration_temperature_tolerance_C': 25, 'filter_aid_preincineration_time_h': 4, 'drying_temperature_C': 130, 'drying_temperature_tolerance_C': 2, 'minimum_drying_time_h': 2, 'incineration_temperature_C': 500, 'incineration_temperature_tolerance_C': 25, 'minimum_incineration_time_h': 2, 'minimum_full_process_blanks_per_batch': 2}.items():
        require(fiber[key] == value, 'Fibre metadata mismatch: ' + key)
    require(fiber['filter_pore_range_um'] == [40, 100], 'Fibre filter pore range')
    require(fiber['repeated_drying_to_constant_mass'] is True and fiber['repeated_incineration_to_constant_mass'] is True, 'Fibre constant mass settings')
    require('mean' in fiber['blank_correction'] and 'blank' in fiber['result_calculation'], 'Fibre full-process blank correction')
    with (ROOT / 'experimental_data/reference_method_metadata_publication_20261014.csv').open(encoding='utf-8', newline='') as stream:
        csv_rows = list(csv.DictReader(stream))
    require(len(csv_rows) == 4, 'Current assay CSV coverage')
    for row in csv_rows:
        current = by_target[row['target']]
        details = json.loads(row['details_json'])
        for key, value in details.items():
            require(current[key] == value, 'Current method CSV/JSON mismatch: ' + key)
        require(row['reporting_basis'] == current['reporting_basis'], 'Current method CSV reporting basis')
    return by_target

def check_acquisition():
    row = read_json('experimental_data/spectral_acquisition_metadata_publication_20261014.json')
    require(row['metadata_version'] == VERSION and row['confirmation_date'] == '2026-10-14', 'Acquisition metadata version')
    for key, value in {'independent_samples': 98, 'ground_sieve_mesh': 40, 'ground_sieve_aperture_mm': 0.425, 'wavelength_start_nm': 400.0, 'wavelength_end_nm': 2499.5, 'wavelength_sampling_interval_nm': 0.5, 'wavelength_variables': 4200, 'internal_subscans_per_measurement': 32, 'independent_repackings_per_sample_per_form': 2, 'exported_spectra_per_sample_per_form': 1, 'modelling_grid_interval_nm': 2.0, 'modelling_grid_end_nm': 2498.0, 'modelling_grid_variables': 1050}.items():
        require(row[key] == value, 'Acquisition metadata mismatch: ' + key)
    require(row['measurement_mode'] == 'diffuse reflectance' and row['exported_spectral_scale'] == 'apparent absorbance', 'Acquisition scale or mode')
    require(row['chemometric_preprocessing_before_export'] == [] and row['raw_export_values_changed'] is False, 'Pre-export treatment or raw values')
    require(row['individual_repacking_spectra_in_archive'] is False, 'Individual spectrum availability')
    require('arithmetic mean' in row['aggregation'] and 'mean-spectrum' in row['legacy_replicate_column_meaning'], 'Exported averaging definition')
    return row

def check_preserved_research():
    preservation = read_json('provenance/Preserved_Research_File_Hashes_20261014_v2.json')
    require(preservation['file_count'] == len(preservation['sha256']), 'Research preservation inventory coverage')
    for relative, expected in preservation['sha256'].items():
        require(release.sha(ROOT / relative) == expected, 'Research file changed: ' + relative)
    correction = read_json('provenance/Publication_Correction_Audit_20261014_v3.json')
    require(correction['raw_data_changed'] is False and correction['frozen_results_changed'] is False and correction['scientific_model_code_changed'] is False, 'Correction scope')
    require(correction['raw_assay_or_spectral_repeat_values_fabricated'] is False and correction['model_retraining'] is False, 'Unsupported correction claim')
    return preservation['file_count']

def check_citation():
    metadata = read_json('PUBLICATION_METADATA.json')
    require(metadata['publication_package_version'] == VERSION and metadata['dataset_title'] == TITLE, 'Package citation metadata')
    require(metadata['status'] == 'local upload package' and metadata['public_upload_performed_by_this_revision'] is False, 'Publication status')
    cff = (ROOT / 'CITATION.cff').read_text(encoding='utf-8')
    require('title: "' + TITLE + '"' in cff and 'version: "' + VERSION + '"' in cff, 'CFF citation mismatch')
    require(cff.count('family-names:') == 8, 'CFF author coverage')
    return metadata

def main():
    coverage = check_inventory_coverage()
    checked = release.check_release_files()
    require(coverage == checked, 'Release checksum inventory count')
    chemistry = release.check_inputs()
    bound = release.check_original_binding()
    predictions, metrics = release.check_predictions(chemistry)
    candidates = release.check_selection()
    selective = release.check_selective()
    methods = check_publication_methods()
    acquisition = check_acquisition()
    preserved = check_preserved_research()
    metadata = check_citation()
    report = {'status': 'passed', 'scope': 'current local publication package', 'version': VERSION, 'dataset_title': TITLE, 'released_files_checked': checked, 'preserved_research_files_checked': preserved, 'samples': len(chemistry), 'reference_constituents': len(methods), 'reference_basis': 'as-is (wet)', 'nitrogen_to_protein_factor': 6.25, 'moisture_comparison_experiment_performed': methods['moisture']['oven_method_validation_completed'], 'independent_repackings_averaged_per_exported_spectrum': acquisition['independent_repackings_per_sample_per_form'], 'internal_subscans_per_measurement': acquisition['internal_subscans_per_measurement'], 'ground_sieve_aperture_mm': acquisition['ground_sieve_aperture_mm'], 'oof_rows': len(predictions), 'recomputed_repeat_metrics': len(metrics), 'inner_candidate_rows': candidates, 'historically_bound_code_and_config_files': bound, 'selective_inferential_comparisons': len(selective), 'raw_and_prepared_numeric_values_match': True, 'nested_split_isolation': True, 'current_author_confirmed_method_metadata': 'passed', 'historical_strict_metadata_scope': 'archived template fields; unchanged', 'raw_laboratory_measurements_included': False, 'individual_repacking_spectra_included': False, 'public_code_comments': 0, 'public_code_docstrings': 0, 'chinese_text_or_filenames': 0, 'external_validation': False, 'full_formal_refit_executed_for_this_release': False, 'public_upload_performed': metadata['public_upload_performed_by_this_revision']}
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
