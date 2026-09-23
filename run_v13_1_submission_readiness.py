#!/usr/bin/env python3
from pathlib import Path
import os, subprocess, sys
for _thread_var in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(_thread_var,'1')
ROOT=Path(__file__).resolve().parent; env=os.environ.copy(); env['PYTHONDONTWRITEBYTECODE']='1'
# Submission readiness is intentionally stricter than computational readiness:
# it requires verified metadata AND a completed, freeze-bound four-target formal run.
subprocess.run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_package.py'),'--root',str(ROOT),'--verify-checksums','--strict-submission-metadata'],check=True,env=env)
subprocess.run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_results.py'),'--root',str(ROOT),'--results-dir',str(ROOT/'v13_2_results'),'--expected-targets','crude_protein','moisture','crude_ash','crude_fiber','--formal'],check=True,env=env)
print('V13.2.0-98 submission-readiness checks PASSED.')
