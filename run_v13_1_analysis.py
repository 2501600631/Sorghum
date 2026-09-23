#!/usr/bin/env python3
"""Immutable formal V13.2.0-98 four-target runner.

This wrapper intentionally accepts no user arguments. It validates frozen data
and active-code hashes before and after analysis, refuses a non-empty formal
results directory, runs the fixed four-target publication design, writes a
freeze receipt binding results to code/data/metadata hashes, and performs strict
result-integrity validation.
"""
from pathlib import Path
import hashlib, json, os, shutil, subprocess, sys
for _thread_var in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(_thread_var,'1')
ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'v13_2_results'
TARGETS=['crude_protein','moisture','crude_ash','crude_fiber']
ENV=os.environ.copy(); ENV['PYTHONDONTWRITEBYTECODE']='1'

def sha256_file(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def run(cmd): subprocess.run(cmd,check=True,env=ENV)

def main():
    if len(sys.argv)!=1:
        raise SystemExit('FORMAL_RUNNER_LOCK: run_v13_1_analysis.py accepts no command-line overrides. Use the smoke/primary-check runners for non-formal work.')
    run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_package.py'),'--root',str(ROOT),'--verify-checksums'])
    if RESULTS.exists() and any(RESULTS.iterdir()):
        raise RuntimeError(f'FORMAL_RESULTS_DIRECTORY_NOT_CLEAN: {RESULTS}. Move/archive/delete prior contents before the one formal run.')
    if RESULTS.exists(): shutil.rmtree(RESULTS)
    run([sys.executable,str(ROOT/'experimental_code'/'sorghum_v13_1_frozen.py'),'--data-dir',str(ROOT/'experimental_data'),'--results-dir',str(RESULTS),'--profile','publication','--targets',*TARGETS])
    run([sys.executable,str(ROOT/'experimental_code'/'v13_1_privileged_utility_analysis.py'),'--results-dir',str(RESULTS),'--profile','publication'])
    run([sys.executable,str(ROOT/'experimental_code'/'v13_1_statistical_robustness.py'),'--data-dir',str(ROOT/'experimental_data'),'--results-dir',str(RESULTS),'--profile','publication'])
    frozen_code=json.loads((ROOT/'FROZEN_CODE_HASHES_V13_1.json').read_text(encoding='utf-8')).get('sha256',{})
    receipt={'framework_version':'13.2.0-98','formal_targets':TARGETS,'data_and_metadata_sha256':{},'active_code_and_config_sha256':{}}
    for rel in ['whole_grain_replicate_spectra.csv.gz','ground_grain_replicate_spectra.csv.gz','reference_chemistry.csv','sample_metadata_template.csv','reference_method_metadata_template.csv','FROZEN_DATA_HASHES.json']:
        p=ROOT/'experimental_data'/rel
        if p.exists(): receipt['data_and_metadata_sha256'][rel]=sha256_file(p)
    for rel in frozen_code:
        p=ROOT/rel
        if p.exists(): receipt['active_code_and_config_sha256'][rel]=sha256_file(p)
    audit=RESULTS/'00_audit'; audit.mkdir(parents=True,exist_ok=True)
    (audit/'formal_run_freeze_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_results.py'),'--root',str(ROOT),'--results-dir',str(RESULTS),'--expected-targets',*TARGETS,'--formal'])
    run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_package.py'),'--root',str(ROOT),'--verify-checksums'])
    print('V13.2.0-98 FORMAL FOUR-TARGET ANALYSIS COMPLETED AND FROZEN:',RESULTS)

if __name__=='__main__': main()
