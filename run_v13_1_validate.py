from pathlib import Path
import os, subprocess, sys
for _thread_var in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(_thread_var,'1')
ROOT=Path(__file__).resolve().parent
env=os.environ.copy(); env['PYTHONDONTWRITEBYTECODE']='1'
subprocess.run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_package.py'),'--root',str(ROOT),*sys.argv[1:]],check=True,env=env)
