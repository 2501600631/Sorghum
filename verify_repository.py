import ast
import gzip
import hashlib
import io
import json
import re
import sys
import tokenize
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'v13_2_results'
TARGETS = ['crude_protein', 'moisture', 'crude_ash', 'crude_fiber']
METHODS = {'nested_strong_whole', 'nested_ground_assisted', 'ground_reference', 'whole_pls', 'whole_ridge', 'whole_elastic', 'whole_ipls', 'whole_svr', 'whole_krr', 'whole_pls_krr', 'whole_extratrees', 'epo_pls', 'epo_elastic', 'pooled_elastic', 'dipls'}
HAN = re.compile('[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]')

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()

def require(condition, message):
    if not condition:
        raise AssertionError(message)

def check_release_files():
    count = 0
    for line in (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        if not line:
            continue
        expected, relative = line.split('  ', 1)
        require(sha(ROOT / relative) == expected, 'Release checksum mismatch: ' + relative)
        count += 1
    for path in ROOT.rglob('*'):
        if not path.is_file() or '__pycache__' in path.parts or '.git' in path.parts:
            continue
        require(not HAN.search(path.relative_to(ROOT).as_posix()), 'Non-English filename: ' + path.name)
        if path.suffix in {'.py', '.json', '.csv', '.txt', '.md', '.bat', '.yml', '.cff', '.svg'} or path.name == '.gitignore':
            text = path.read_text(encoding='utf-8-sig')
            require(not HAN.search(text), 'Chinese text remains: ' + str(path.relative_to(ROOT)))
            require(not re.search(r'(?:[A-Z]:[\\/]+Users[\\/]+\d{3,}[\\/]|api_key\s*[:=]|password\s*[:=])', text, re.I), 'Potential local path or credential: ' + path.name)
            if path.suffix == '.py':
                require(not any(t.type == tokenize.COMMENT for t in tokenize.generate_tokens(io.StringIO(text).readline)), 'Code comment remains: ' + path.name)
                tree = ast.parse(text)
                for node in ast.walk(tree):
                    if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                        require(ast.get_docstring(node) is None, 'Code docstring remains: ' + path.name)
                compile(text, str(path), 'exec')
        elif path.name.endswith('.csv.gz'):
            with gzip.open(path, 'rt', encoding='utf-8-sig') as stream:
                require(not HAN.search(stream.read()), 'Chinese text remains in compressed CSV: ' + path.name)
        elif path.suffix == '.xlsx':
            with zipfile.ZipFile(path) as book:
                for entry in book.namelist():
                    if entry.endswith('.xml'):
                        require(not HAN.search(book.read(entry).decode('utf-8')), 'Chinese workbook metadata: ' + path.name)
    return count

def check_inputs():
    chemistry = pd.read_csv(ROOT / 'experimental_data/reference_chemistry.csv')
    require(len(chemistry) == 98 and chemistry.sample_id.nunique() == 98, 'Reference sample count')
    require(np.isfinite(chemistry[TARGETS].to_numpy(float)).all(), 'Nonfinite chemistry')
    raw = pd.read_excel(ROOT / 'raw_inputs/reference_chemistry_raw.xlsx')
    require(raw.sample_id.astype(str).tolist() == chemistry.sample_id.astype(str).tolist(), 'Chemistry IDs or row order changed')
    raw_columns = ['moisture_pct', 'crude_protein_pct', 'crude_ash_pct', 'crude_fiber_pct']
    np.testing.assert_allclose(raw[raw_columns].to_numpy(float), chemistry[['moisture', 'crude_protein', 'crude_ash', 'crude_fiber']].to_numpy(float), rtol=0, atol=1e-12)
    for form in ['whole', 'ground']:
        prepared = pd.read_csv(ROOT / 'experimental_data' / (form + '_grain_replicate_spectra.csv.gz'))
        sheet = form + '_grain_400_2499p5'
        spectral = pd.read_excel(ROOT / 'raw_inputs/paired_spectra_raw.xlsx', sheet_name=sheet)
        require(prepared.shape == (98, 4202) and spectral.shape == (98, 4201), 'Spectral dimensions')
        require(prepared.replicate.eq(1).all(), 'Unexpected technical replicates')
        require(prepared.sample_id.astype(str).tolist() == spectral.sample_id.astype(str).tolist(), 'Spectrum IDs or row order changed')
        require(set(prepared.sample_id.astype(str)) == set(chemistry.sample_id.astype(str)), 'Paired input ID mismatch')
        np.testing.assert_array_equal(np.asarray(spectral.columns[1:], dtype=float), np.arange(400, 2500, 0.5))
        np.testing.assert_allclose(prepared.iloc[:, 2:].to_numpy(float), spectral.iloc[:, 1:].to_numpy(float), rtol=0, atol=1e-12)
        require(np.isfinite(prepared.iloc[:, 2:].to_numpy(float)).all(), 'Nonfinite spectral values')
    require(float(chemistry.set_index('sample_id').loc['202512-06142', 'crude_protein']) == 9.2, 'Frozen protein record changed')
    return chemistry

def check_original_binding():
    mapping = json.loads((ROOT / 'provenance/Publication_Transformation.json').read_text(encoding='utf-8'))
    by_original = {row['original_path']: row for row in mapping['files']}
    receipt = json.loads((RESULTS / '00_audit/formal_run_freeze_receipt.json').read_text(encoding='utf-8'))
    for relative, historical_hash in receipt['active_code_and_config_sha256'].items():
        row = by_original[relative]
        require(row['original_sha256'] == historical_hash, 'Original code binding mismatch: ' + relative)
        require(sha(ROOT / row['public_path']) == row['public_sha256'], 'Public derivative mismatch: ' + relative)
        if relative.endswith('.py'):
            require(row.get('scientific_executable_ast_preserved') is True, 'Executable AST mismatch: ' + relative)
    for relative, historical_hash in receipt['data_and_metadata_sha256'].items():
        path = ROOT / 'experimental_data' / relative
        if relative == 'FROZEN_DATA_HASHES.json':
            path = ROOT / 'provenance/Original_Data_Freeze_Hashes.json'
        require(sha(path) == historical_hash, 'Original input binding mismatch: ' + relative)
    for row in mapping['files']:
        if row['original_path'].startswith('v13_2_results/'):
            require(sha(ROOT / row['public_path']) == row['original_sha256'], 'Frozen result changed: ' + row['public_path'])
    return len(receipt['active_code_and_config_sha256'])

def check_predictions(chemistry):
    oof = pd.read_csv(RESULTS / '01_predictions/all_targets_outer_oof_predictions.csv')
    require(len(oof) == 58800 and set(oof.target) == set(TARGETS) and set(oof.method) == METHODS, 'OOF dimensions or model coverage')
    require(oof.sample_id.nunique() == 98 and oof.outer_repeat.nunique() == 10 and oof.outer_fold.nunique() == 5, 'OOF sample/repeat/fold count')
    require(not oof.duplicated(['target', 'method', 'sample_id', 'outer_repeat']).any(), 'Duplicate OOF records')
    require(oof.groupby(['target', 'method', 'sample_id']).size().eq(10).all(), 'Incomplete OOF repetitions')
    require(oof.groupby(['target', 'method', 'outer_repeat']).sample_id.nunique().eq(98).all(), 'Incomplete repeat coverage')
    require(np.isfinite(oof[['y_true', 'y_pred']].to_numpy(float)).all(), 'Nonfinite OOF values')
    reference = chemistry.set_index('sample_id')
    for target, rows in oof.groupby('target'):
        np.testing.assert_allclose(reference.loc[rows.sample_id, target].to_numpy(float), rows.y_true.to_numpy(float), rtol=0, atol=1e-12)
    saved = pd.read_csv(RESULTS / '03_performance/outer_repeat_performance.csv')
    metrics = []
    for keys, rows in oof.groupby(['target', 'method', 'outer_repeat']):
        true = rows.y_true.to_numpy(float)
        predicted = rows.y_pred.to_numpy(float)
        squared = np.sum((true - predicted) ** 2)
        metrics.append({'target': keys[0], 'method': keys[1], 'outer_repeat': keys[2], 'rmse_calc': np.sqrt(squared / len(rows)), 'r2_calc': 1 - squared / np.sum((true - true.mean()) ** 2)})
    merged = pd.DataFrame(metrics).merge(saved, on=['target', 'method', 'outer_repeat'], validate='one_to_one')
    require(len(merged) == 600, 'Saved repeat metric count')
    np.testing.assert_allclose(merged.rmse_calc, merged.rmse, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(merged.r2_calc, merged.r2, rtol=1e-12, atol=1e-12)
    return oof, merged

def check_selection():
    selected = pd.read_csv(RESULTS / '02_selection/outer_fold_selected_workflows.csv')
    require(len(selected) == 200, 'Selected workflow count')
    for row in selected.itertuples():
        train = set(json.loads(row.outer_train_sample_ids_json))
        test = set(json.loads(row.outer_test_sample_ids_json))
        require(not train & test and len(train | test) == 98, 'Outer-fold overlap or incomplete coverage')
        for split in json.loads(row.inner_split_sample_ids_json):
            inner_train = set(split['train_sample_ids'])
            inner_test = set(split['validation_sample_ids'])
            require(not inner_train & inner_test and inner_train | inner_test == train and not (inner_train | inner_test) & test, 'Inner-fold overlap or outer-test leakage')
    require(selected.representation_candidates_per_branch.eq(30).all() and selected.representations_retained_per_branch.eq(2).all(), 'Representation count changed')
    screen = pd.read_csv(RESULTS / '02_selection/representation_screening_trace.csv')
    require(len(screen) == 18000 and screen.screen_candidate_count.eq(4).all() and screen.all_screen_proxies_valid.astype(bool).all(), 'Screening coverage')
    candidates = pd.read_csv(RESULTS / '02_selection/inner_cv_candidate_trace.csv.gz')
    require(len(candidates) == 152400 and candidates.fit_status.eq('valid').all(), 'Candidate coverage or failure status')
    return len(candidates)

def check_selective():
    inference = pd.read_csv(RESULTS / '06_robustness/selective_grinding_vs_random_inference.csv')
    require(len(inference) == 8 and inference.selection_scope.eq('within_outer_fold').all(), 'Selective comparison count or scope')
    require(set(np.round(inference.requested_budget, 8)) == {0.1, 0.2}, 'Selective inferential budgets')
    expected = {'crude_protein': (0.20545205639016575, 0.22226175491402564, 0.04409779511024449), 'crude_fiber': (0.20318197987073483, 0.24682267170797406, 0.0007999600019999)}
    rows = inference[inference.requested_budget.eq(0.2)].set_index('target')
    require(set(rows.index) == set(TARGETS), 'Missing 20-percent comparisons')
    for target, values in expected.items():
        actual = rows.loc[target, ['guided_rmse', 'matched_random_expected_rmse', 'holm_adjusted_p_all_target_budget_tests']].to_numpy(float)
        np.testing.assert_allclose(actual, np.asarray(values), rtol=0, atol=1e-12)
    return inference

def main():
    released_files = check_release_files()
    chemistry = check_inputs()
    bound_files = check_original_binding()
    predictions, metrics = check_predictions(chemistry)
    candidates = check_selection()
    selective = check_selective()
    report = {'status': 'passed', 'released_files_checked': released_files, 'samples': len(chemistry), 'oof_rows': len(predictions), 'recomputed_repeat_metrics': len(metrics), 'inner_candidate_rows': candidates, 'historically_bound_code_and_config_files': bound_files, 'selective_inferential_comparisons': len(selective), 'raw_and_prepared_numeric_values_match': True, 'nested_split_isolation': True, 'public_code_comments': 0, 'public_code_docstrings': 0, 'chinese_text_or_filenames': 0, 'external_validation': False, 'full_formal_refit_executed_for_this_release': False}
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    sys.dont_write_bytecode = True
    main()
