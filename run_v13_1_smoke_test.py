#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, sys
for _thread_var in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(_thread_var,'1')
ROOT=Path(__file__).resolve().parent; out=ROOT/'v13_2_smoke_results'; env=os.environ.copy(); env['PYTHONDONTWRITEBYTECODE']='1'
if out.exists(): shutil.rmtree(out)
subprocess.run([sys.executable,str(ROOT/'experimental_code'/'sorghum_v13_1_frozen.py'),'--data-dir',str(ROOT/'experimental_data'),'--results-dir',str(out),'--profile','smoke','--targets','crude_protein','moisture','crude_ash','crude_fiber','--outer-splits','2','--outer-repeats','1','--inner-splits','2','--bootstrap-repeats','20','--signflip-repeats','20','--top-representations-per-branch','1','--smoke-sample-limit','98'],check=True,env=env)
subprocess.run([sys.executable,str(ROOT/'experimental_code'/'v13_1_privileged_utility_analysis.py'),'--results-dir',str(out),'--profile','smoke','--seed','20260830','--random-repeats','5','--bootstrap-repeats','20'],check=True,env=env)
subprocess.run([sys.executable,str(ROOT/'experimental_code'/'v13_1_statistical_robustness.py'),'--data-dir',str(ROOT/'experimental_data'),'--results-dir',str(out),'--profile','smoke','--bootstrap-repeats','20','--signflip-repeats','50','--label-permutations','200'],check=True,env=env)
subprocess.run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_results.py'),'--root',str(ROOT),'--results-dir',str(out),'--expected-targets','crude_protein','moisture','crude_ash','crude_fiber'],check=True,env=env)
print('V13.2.0-98 end-to-end smoke test complete and result integrity passed:',out)
