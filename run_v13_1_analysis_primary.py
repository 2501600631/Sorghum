from pathlib import Path
import os, shutil, subprocess, sys
for _thread_var in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ.setdefault(_thread_var,'1')
ROOT=Path(__file__).resolve().parent; OUT=ROOT/'v13_1_primary_check_results'; ENV=os.environ.copy(); ENV['PYTHONDONTWRITEBYTECODE']='1'
def run(cmd): subprocess.run(cmd,check=True,env=ENV)
def main():
    if len(sys.argv)!=1: raise SystemExit('Primary-check runner accepts no overrides.')
    run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_package.py'),'--root',str(ROOT),'--verify-checksums'])
    if OUT.exists(): shutil.rmtree(OUT)
    run([sys.executable,str(ROOT/'experimental_code'/'sorghum_v13_1_frozen.py'),'--data-dir',str(ROOT/'experimental_data'),'--results-dir',str(OUT),'--profile','publication','--targets','crude_protein','moisture'])
    run([sys.executable,str(ROOT/'experimental_code'/'v13_1_privileged_utility_analysis.py'),'--results-dir',str(OUT),'--profile','publication'])
    run([sys.executable,str(ROOT/'experimental_code'/'v13_1_statistical_robustness.py'),'--data-dir',str(ROOT/'experimental_data'),'--results-dir',str(OUT),'--profile','publication'])
    run([sys.executable,str(ROOT/'experimental_code'/'validate_v13_1_results.py'),'--root',str(ROOT),'--results-dir',str(OUT),'--expected-targets','crude_protein','moisture'])
    (OUT/'NOT_FOR_SUBMISSION.txt').write_text('Two-target preflight only. Formal manuscript results must come from run_v13_1_analysis.py with all four locked targets.\n',encoding='utf-8')
    print('V13.2.0-98 PRIMARY CHECK COMPLETED (NOT FOR SUBMISSION):',OUT)
if __name__=='__main__': main()
