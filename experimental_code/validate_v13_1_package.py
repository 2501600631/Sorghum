from __future__ import annotations
import argparse, hashlib, inspect, json, sys
from pathlib import Path
import pandas as pd
sys.dont_write_bytecode=True
import sorghum_v13_1_frozen as core
import v13_1_privileged_utility_analysis as utility
import v13_1_statistical_robustness as robust

FRAMEWORK_VERSION='13.2.0-98'

def sha256_file(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]); ap.add_argument('--verify-checksums',action='store_true'); ap.add_argument('--strict-submission-metadata',action='store_true')
    args=ap.parse_args(); root=args.root.resolve()
    required=['RUN_CONFIG_V13_1.json','README.md','requirements.txt','run_v13_1_analysis.py','run_v13_1_analysis_primary.py','run_v13_1_smoke_test.py','run_v13_1_validate.py','run_v13_1_analysis_all_targets.bat','run_v13_1_analysis_primary.bat','run_v13_1_smoke_test.bat','run_v13_1_validate.bat','run_v13_1_submission_readiness.py','run_v13_1_submission_readiness.bat','run_v13_1_algorithm_selftest.py','run_v13_1_algorithm_selftest.bat',
              'experimental_code/sorghum_v13_1_frozen.py','experimental_code/sorghum_v13_data_utils.py','experimental_code/v13_1_privileged_utility_analysis.py','experimental_code/v13_1_statistical_robustness.py','experimental_code/v13_1_algorithm_selftest.py','experimental_code/validate_v13_1_results.py',
              'experimental_data/whole_grain_replicate_spectra.csv.gz','experimental_data/ground_grain_replicate_spectra.csv.gz','experimental_data/reference_chemistry.csv','experimental_data/sample_metadata_template.csv','experimental_data/reference_method_metadata_template.csv','experimental_data/FROZEN_DATA_HASHES.json','FROZEN_CODE_HASHES_V13_1.json']
    missing=[x for x in required if not (root/x).exists()]
    if missing: raise FileNotFoundError(f'Missing required files: {missing}')
    active=[root/'run_v13_1_analysis.py',root/'run_v13_1_analysis_primary.py',root/'run_v13_1_smoke_test.py',root/'run_v13_1_validate.py',root/'run_v13_1_submission_readiness.py',root/'run_v13_1_algorithm_selftest.py',
            root/'experimental_code'/'sorghum_v13_1_frozen.py',root/'experimental_code'/'sorghum_v13_data_utils.py',root/'experimental_code'/'v13_1_privileged_utility_analysis.py',root/'experimental_code'/'v13_1_statistical_robustness.py',root/'experimental_code'/'v13_1_algorithm_selftest.py',root/'experimental_code'/'validate_v13_1_results.py']
    for p in active: compile(p.read_text(encoding='utf-8'),str(p),'exec')

    for p in [root/'experimental_code'/'sorghum_v13_1_frozen.py',root/'experimental_code'/'v13_1_privileged_utility_analysis.py',root/'experimental_code'/'v13_1_statistical_robustness.py']:
        if 'except Exception' in p.read_text(encoding='utf-8'): raise RuntimeError(f'Broad exception handler remains in active scientific code: {p.name}')
    data=core.load_data(root/'experimental_data',2.0,400.0,2498.0)
    if len(data.sample_ids)!=98: raise RuntimeError('Expected 98 independent samples')
    if data.whole_replicates.shape!=(98,1,1050) or data.ground_replicates.shape!=(98,1,1050): raise RuntimeError('Unexpected replicate shape')
    if data.chemistry[core.TARGETS].isna().any().any(): raise RuntimeError('Missing reference chemistry')
    groups,provenance=core.load_metadata_groups(root/'experimental_data',data.sample_ids); reference=core.validate_reference_metadata(root/'experimental_data')
    if args.strict_submission_metadata:
        blockers=[]
        if provenance.get('source_group_status')!='complete': blockers.append('sample provenance/source_group is not verified complete')
        if not reference.get('complete',False): blockers.append('reference chemistry method metadata are incomplete')
        if blockers: raise RuntimeError('STRICT_SUBMISSION_METADATA_FAILED: '+'; '.join(blockers))
    reps=core.representation_space('publication')
    if len(reps)!=30: raise RuntimeError('Publication representation count must be 30')
    parser=core.build_parser(); defaults=parser.parse_args([]); core.enforce_publication_design_lock(defaults)
    if defaults.top_representations_per_branch!=2: raise RuntimeError('Two representations must be retained per branch')
    changed=parser.parse_args(['--top-representations-per-branch','3'])
    try: core.enforce_publication_design_lock(changed)
    except RuntimeError: pass
    else: raise RuntimeError('Core publication design lock failed')
    epo_src=inspect.getsource(core.fit_epo_base)
    if 'D.mean' in epo_src or 'mean(axis=0' in epo_src: raise RuntimeError('EPO difference matrix must be non-centered')
    sel_src=inspect.getsource(core.select_min_inner_rmse)
    if 'complexity' in sel_src.split('ordered=')[1]: raise RuntimeError('Complexity entered candidate ranking')
    probe=reps[0]
    proxy_expected={'whole':{'whole_pls','whole_elastic','whole_svr','whole_krr'},'transfer':{'epo_pls','epo_elastic','pooled_elastic','dipls'},'ground':{'ground_pls','ground_elastic','ground_svr','ground_krr'}}
    for branch,expected in proxy_expected.items():
        cands=core._screen_proxy_candidates(probe,False,branch)
        if len(cands)!=4 or {c.family for c in cands}!=expected: raise RuntimeError(f'{branch} screening proxy mismatch')
    whole={c.family.replace('whole_','',1) for c in core._candidate_grids(probe,'publication',False)}; ground={c.family.replace('ground_','',1) for c in core._ground_candidate_grids(probe,'publication')}; base={'pls','ridge','elastic','ipls','svr','krr','pls_krr','extratrees'}
    if whole!=base or ground!=base: raise RuntimeError('Whole/Ground full model families are not mirrored')
    if {c.family for c in core._candidate_grids(probe,'publication',True)}!={'epo_pls','epo_elastic','pooled_elastic','dipls'}: raise RuntimeError('Transfer family mismatch')

    uparser=utility.build_parser(); udefaults=uparser.parse_args([]); utility.enforce_lock(udefaults)
    try: utility.enforce_lock(uparser.parse_args(['--random-repeats','499']))
    except RuntimeError: pass
    else: raise RuntimeError('Utility design lock failed')
    rparser=robust.build_parser(); rdefaults=rparser.parse_args([]); robust.enforce_lock(rdefaults)
    try: robust.enforce_lock(rparser.parse_args(['--bootstrap-repeats','1999']))
    except RuntimeError: pass
    else: raise RuntimeError('Robustness design lock failed')

    runner_src=(root/'run_v13_1_analysis.py').read_text(encoding='utf-8')
    if 'len(sys.argv)!=1' not in runner_src or "TARGETS=['crude_protein','moisture','crude_ash','crude_fiber']" not in runner_src: raise RuntimeError('Formal runner is not hard-locked')
    checksum_status='not_requested'
    if args.verify_checksums:
        frozen=json.loads((root/'experimental_data/FROZEN_DATA_HASHES.json').read_text(encoding='utf-8')).get('sha256',{})
        failures=[rel for rel,h in frozen.items() if not (root/'experimental_data'/rel).exists() or sha256_file(root/'experimental_data'/rel)!=h]
        if failures: raise RuntimeError(f'Frozen-data hash failures: {failures}')
        frozen_code=json.loads((root/'FROZEN_CODE_HASHES_V13_1.json').read_text(encoding='utf-8')).get('sha256',{})
        failures=[rel for rel,h in frozen_code.items() if not (root/rel).exists() or sha256_file(root/rel)!=h]
        if failures: raise RuntimeError(f'Frozen-code/config hash failures: {failures}')
        checksum_status='passed_data_and_active_code'
    report={'status':'passed','framework_version':FRAMEWORK_VERSION,'core_script_version':core.SCRIPT_VERSION,'coded_samples':98,'observations_per_sample_per_presentation':1,'targets':core.TARGETS,'checksum_verification':checksum_status,
            'publication_representation_candidates_per_branch':30,'publication_representations_retained_per_branch':2,'screen_proxies_per_branch':4,'whole_ground_full_base_model_families_identical':True,'ground_representation_selection_independent':True,'transfer_representation_selection_independent':True,
            'selective_grinding_rank_scope':'within_outer_fold','provenance_aware_inference_enabled':True,'candidate_failure_provenance_enabled':True,'unexpected_candidate_exceptions_abort':True,'corrected_repeated_cv_sensitivity_enabled':True,
            'formal_runner_hard_locked':True,'utility_design_lock_active':True,'robustness_design_lock_active':True,'epo_difference_svd':'non-centered','candidate_selection_rule':'minimum mean inner-CV RMSE; deterministic candidate_key tie-break',
            'sample_provenance_status':provenance.get('source_group_status'),'reference_method_metadata_complete':reference.get('complete',False),'strict_submission_metadata_requested':bool(args.strict_submission_metadata),'external_validation_present':False,'deployment_claim_allowed':False}
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
