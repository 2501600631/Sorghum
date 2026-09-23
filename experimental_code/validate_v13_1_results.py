#!/usr/bin/env python3
"""Strict completeness and freeze-binding validator for V13.2.0-98 results."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

FRAMEWORK_VERSION='13.2.0-98'
ALL_TARGETS=['crude_protein','moisture','crude_ash','crude_fiber']
EXPECTED_METHODS={'nested_strong_whole','nested_ground_assisted','ground_reference','whole_pls','whole_ridge','whole_elastic','whole_ipls','whole_svr','whole_krr','whole_pls_krr','whole_extratrees','epo_pls','epo_elastic','pooled_elastic','dipls'}
ALLOWED_INVALID_FAILURES={'ElasticNetConvergenceError','LinAlgError','FloatingPointError','NonFinitePrediction','IncompleteInnerCV'}

def require(p:Path)->Path:
    if not p.exists(): raise FileNotFoundError(p)
    return p

def sha256_file(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]); ap.add_argument('--results-dir',type=Path,required=True); ap.add_argument('--expected-targets',nargs='+',default=ALL_TARGETS); ap.add_argument('--formal',action='store_true')
    args=ap.parse_args(); root=args.root.resolve(); r=args.results_dir.resolve(); targets=list(args.expected_targets)
    if args.formal and targets!=ALL_TARGETS: raise RuntimeError(f'Formal V13.2.0-98 targets are locked to {ALL_TARGETS}')
    files=['ANALYSIS_SUMMARY.json','00_audit/run_manifest.json','00_audit/sample_provenance_audit.json','00_audit/model_space_fairness_audit.json','00_audit/numerical_convergence_summary.json',
           '01_predictions/all_targets_outer_oof_predictions.csv','02_selection/outer_fold_selected_workflows.csv','02_selection/representation_screening_trace.csv','02_selection/inner_cv_candidate_trace.csv.gz','02_selection/outer_test_split_manifest.csv',
           '03_performance/method_performance_summary.csv','04_inference/ground_assisted_vs_strong_whole_cluster_bootstrap.csv','04_inference/ground_assisted_vs_strong_whole_signflip.csv',
           '05_utility_boundary/trait_privileged_incremental_value.csv','05_utility_boundary/selective_grinding_cost_performance_summary.csv','05_utility_boundary/selective_grinding_sample_decisions.csv.gz','05_utility_boundary/uncertainty_validity_summary.csv',
           '06_robustness/technical_replicate_status.json','06_robustness/practical_gain_sensitivity.csv','06_robustness/selective_grinding_vs_random_inference.csv','06_robustness/corrected_repeated_cv_sensitivity.csv',
           '06_robustness/primary_protein_oof_label_permutation_summary.json','06_robustness/primary_protein_oof_label_permutation_null.csv.gz']
    if args.formal: files.append('00_audit/formal_run_freeze_receipt.json')
    for f in files: require(r/f)
    summary=json.loads((r/'ANALYSIS_SUMMARY.json').read_text()); manifest=json.loads((r/'00_audit/run_manifest.json').read_text())
    if manifest.get('framework_version')!=FRAMEWORK_VERSION: raise RuntimeError('Run-manifest framework version mismatch')
    if set(summary.get('targets',[]))!=set(targets): raise RuntimeError('Analysis-summary targets mismatch')
    oof=pd.read_csv(r/'01_predictions/all_targets_outer_oof_predictions.csv')
    if set(oof.target.unique())!=set(targets): raise RuntimeError('OOF target mismatch')
    required_cols={'target','method','sample_id','inference_group','outer_repeat','outer_fold','y_true','y_pred'}
    if not required_cols.issubset(oof.columns): raise RuntimeError(f'OOF missing columns: {sorted(required_cols-set(oof.columns))}')
    if oof[['target','method','sample_id','outer_repeat']].duplicated().any(): raise RuntimeError('Duplicate OOF target-method-sample-repeat rows')
    if oof.inference_group.astype(str).str.strip().eq('').any(): raise RuntimeError('Blank inference_group in OOF')
    if not np.isfinite(oof[['y_true','y_pred']].to_numpy(float)).all(): raise RuntimeError('Non-finite OOF values')
    n_samples=int(oof.sample_id.nunique()); n_repeats=int(oof.outer_repeat.nunique())
    if args.formal:
        if n_samples!=98 or n_repeats!=10 or int(oof.outer_fold.nunique())!=5: raise RuntimeError('Formal sample/repeat/fold mismatch')
        if set(oof.method.unique())!=EXPECTED_METHODS: raise RuntimeError(f'Formal method set mismatch: got={sorted(oof.method.unique())}')
        expected_rows=4*len(EXPECTED_METHODS)*98*10
        if len(oof)!=expected_rows: raise RuntimeError(f'Formal OOF row count {len(oof)} != {expected_rows}')
        if not (oof.groupby(['target','method','sample_id']).size()==10).all(): raise RuntimeError('Not every target-method-sample has 10 OOF predictions')
    chemistry=pd.read_csv(root/'experimental_data/reference_chemistry.csv').set_index('sample_id')
    expected_ids=set(chemistry.index.astype(str))
    if not set(oof.sample_id.astype(str)).issubset(expected_ids): raise RuntimeError('Unknown OOF sample IDs')
    if args.formal and set(oof.sample_id.astype(str))!=expected_ids: raise RuntimeError('Incomplete formal sample IDs')
    for target,g in oof.groupby('target'):
        actual=chemistry.loc[g.sample_id,target].to_numpy(float)
        if not np.allclose(actual,g.y_true.to_numpy(float),rtol=0,atol=1e-10): raise RuntimeError('OOF reference values changed')
    coverage=oof.groupby(['target','method','outer_repeat']).sample_id.nunique()
    if not (coverage==n_samples).all(): raise RuntimeError('Incomplete OOF coverage')
    sel=pd.read_csv(r/'02_selection/outer_fold_selected_workflows.csv')
    for _,row in sel.iterrows():
        tr=set(json.loads(row.outer_train_sample_ids_json)); te=set(json.loads(row.outer_test_sample_ids_json))
        if tr & te or len(tr|te)!=n_samples: raise RuntimeError('Outer split leakage/coverage failure')
        for inner in json.loads(row.inner_split_sample_ids_json):
            a=set(inner['train_sample_ids']); b=set(inner['validation_sample_ids'])
            if a & b or a|b != tr or (a|b)&te: raise RuntimeError('Inner split leakage/coverage failure')
    if args.formal and len(sel)!=4*50: raise RuntimeError(f'Formal selection rows {len(sel)} != 200')
    if args.formal and (not (sel.representation_candidates_per_branch==30).all() or not (sel.representations_retained_per_branch==2).all()): raise RuntimeError('Formal representation design mismatch')
    screen=pd.read_csv(r/'02_selection/representation_screening_trace.csv')
    if args.formal:
        if len(screen)!=4*50*3*30: raise RuntimeError('Formal screening row count mismatch')
        if not (screen.screen_candidate_count==4).all() or not screen.all_screen_proxies_valid.astype(bool).all(): raise RuntimeError('Formal screening proxy validity failure')
        if set(screen.branch.unique())!={'whole','transfer','ground'}: raise RuntimeError('Missing screening branch')
    tune=pd.read_csv(r/'02_selection/inner_cv_candidate_trace.csv.gz')
    if not {'fit_status','failure_type','failure_message','candidate_key','candidate_group','target','outer_repeat','outer_fold'}.issubset(tune.columns): raise RuntimeError('Inner tuning trace lacks failure provenance')
    if args.formal:
        expected_tune=4*50*2*(126+129+126)
        if len(tune)!=expected_tune: raise RuntimeError(f'Formal inner-tuning rows {len(tune)} != {expected_tune}')
    invalid=tune[tune.fit_status!='valid']
    if len(invalid):
        bad=set(invalid.failure_type.astype(str))-ALLOWED_INVALID_FAILURES
        if bad: raise RuntimeError(f'Unexpected candidate failure types: {sorted(bad)}')
        if invalid.failure_message.fillna('').astype(str).str.len().eq(0).any(): raise RuntimeError('Invalid candidate without failure message')
    # Every globally selected candidate must be valid in the corresponding tuning trace.
    for prefix,group in [('whole','whole'),('transfer','transfer'),('ground','ground_reference')]:
        keycol=f'{prefix}_selected_candidate_key'
        if keycol not in sel.columns: raise RuntimeError(f'Missing {keycol}')
        tmp=sel[['target','outer_repeat','outer_fold',keycol]].merge(tune[tune.candidate_group==group][['target','outer_repeat','outer_fold','candidate_key','fit_status']],left_on=['target','outer_repeat','outer_fold',keycol],right_on=['target','outer_repeat','outer_fold','candidate_key'],how='left',validate='one_to_one')
        if tmp.fit_status.isna().any() or not (tmp.fit_status=='valid').all(): raise RuntimeError(f'Selected {prefix} candidate is missing/invalid in tuning trace')
    numerical=json.loads((r/'00_audit/numerical_convergence_summary.json').read_text())
    if numerical.get('unresolved_selected_final_fits')!=0: raise RuntimeError('Unresolved selected/family ElasticNet final fits')
    fairness=json.loads((r/'00_audit/model_space_fairness_audit.json').read_text())
    if fairness.get('framework_version')!=FRAMEWORK_VERSION or fairness.get('screen_proxy_count_per_branch')!=4: raise RuntimeError('Fairness audit mismatch')
    rs=json.loads((r/'06_robustness/technical_replicate_status.json').read_text())
    if rs.get('status')!='not_applicable' or rs.get('observations_per_sample_per_presentation')!=1:
        raise RuntimeError('Independent-sample replicate-status mismatch')
    if (r/'06_robustness/technical_replicate_oof_predictions.csv').exists():
        raise RuntimeError('Stale technical replicate results in independent-sample run')
    practical=pd.read_csv(r/'06_robustness/practical_gain_sensitivity.csv')
    if len(practical)!=len(targets)*3 or set(np.round(practical.sensitivity_margin_pct,6))!={2.5,5.0,10.0}: raise RuntimeError('Practical-gain output mismatch')
    selective=pd.read_csv(r/'06_robustness/selective_grinding_vs_random_inference.csv')
    if len(selective)!=len(targets)*2 or set(np.round(selective.requested_budget,6))!={0.1,0.2} or not (selective.selection_scope=='within_outer_fold').all(): raise RuntimeError('Selective-grinding inferential output mismatch')
    corrected=pd.read_csv(r/'06_robustness/corrected_repeated_cv_sensitivity.csv')
    if len(corrected)!=len(targets) or (args.formal and not (corrected.n_outer_fold_differences==50).all()): raise RuntimeError('Corrected repeated-CV output mismatch')
    perm=pd.read_csv(r/'06_robustness/primary_protein_oof_label_permutation_null.csv.gz'); ps=json.loads((r/'06_robustness/primary_protein_oof_label_permutation_summary.json').read_text())
    if ps.get('status')=='completed':
        expected_perm=20000 if args.formal else int(ps.get('permutations',len(perm)))
        if len(perm)!=expected_perm: raise RuntimeError('Primary label-permutation count mismatch')
    elif ps.get('status')=='skipped_grouped_provenance':
        if len(perm)!=0: raise RuntimeError('Grouped-provenance permutation should have empty null distribution')
    else: raise RuntimeError('Unknown primary permutation status')
    if args.formal:
        receipt=json.loads((r/'00_audit/formal_run_freeze_receipt.json').read_text())
        if receipt.get('framework_version')!=FRAMEWORK_VERSION or receipt.get('formal_targets')!=ALL_TARGETS: raise RuntimeError('Formal freeze receipt mismatch')
        for rel,expected in receipt.get('data_and_metadata_sha256',{}).items():
            p=root/'experimental_data'/rel
            if not p.exists() or sha256_file(p)!=expected: raise RuntimeError(f'Formal result invalidated by changed data/metadata: {rel}')
        for rel,expected in receipt.get('active_code_and_config_sha256',{}).items():
            p=root/rel
            if not p.exists() or sha256_file(p)!=expected: raise RuntimeError(f'Formal result invalidated by changed active code/config: {rel}')
        frozen_data=json.loads((root/'experimental_data/FROZEN_DATA_HASHES.json').read_text()).get('sha256',{})
        for rel,expected in frozen_data.items():
            if sha256_file(root/'experimental_data'/rel)!=expected: raise RuntimeError(f'Frozen data hash mismatch: {rel}')
        frozen_code=json.loads((root/'FROZEN_CODE_HASHES_V13_1.json').read_text()).get('sha256',{})
        for rel,expected in frozen_code.items():
            if sha256_file(root/rel)!=expected: raise RuntimeError(f'Frozen code/config hash mismatch: {rel}')
        for rel,expected in manifest.get('input_sha256',{}).items():
            if sha256_file(root/'experimental_data'/rel)!=expected: raise RuntimeError(f'Run-manifest input hash mismatch: {rel}')
    report={'status':'passed','framework_version':FRAMEWORK_VERSION,'formal':bool(args.formal),'targets':targets,'samples':n_samples,'outer_repeats':n_repeats,'oof_rows':len(oof),'selection_rows':len(sel),'screening_rows':len(screen),'inner_tuning_rows':len(tune),'invalid_expected_numeric_candidates':len(invalid),'technical_replicate_analysis':'not_applicable','practical_gain_rows':len(practical),'selective_inference_rows':len(selective),'corrected_repeated_cv_rows':len(corrected),'label_permutation_status':ps.get('status')}
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
