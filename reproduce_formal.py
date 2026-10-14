import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Reproduce the fixed formal computation without overwriting historical results.')
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    destination = args.output_dir.resolve()
    if destination.exists():
        raise FileExistsError('The output directory must not exist: ' + str(destination))
    environment = os.environ.copy()
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    with tempfile.TemporaryDirectory(prefix='nir98_formal_') as temporary:
        work = Path(temporary) / 'experiment'
        shutil.copytree(root, work, ignore=shutil.ignore_patterns('v13_2_results', 'v13_2_smoke_results', '__pycache__', '.git', 'rerun_results'))
        subprocess.run([sys.executable, str(work / 'freeze_package.py')], check=True, env=environment)
        subprocess.run([sys.executable, str(work / 'run_v13_1_analysis.py')], check=True, env=environment)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(work / 'v13_2_results', destination)
    print('Completed a new formal computation at ' + str(destination))

if __name__ == '__main__':
    main()
