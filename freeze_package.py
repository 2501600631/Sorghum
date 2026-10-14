from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    results=ROOT/'v13_2_results'
    if results.exists() and any(results.iterdir()):
        raise RuntimeError('Archive the complete formal results directory before refreezing changed inputs/code.')
    d=ROOT/'experimental_data'
    files={p.name:sha(p) for p in sorted(d.iterdir()) if p.is_file() and p.name!='FROZEN_DATA_HASHES.json'}
    (d/'FROZEN_DATA_HASHES.json').write_text(json.dumps({'framework_version':'13.2.0-98','sha256':files},indent=2),encoding='utf-8')
    active=[p for p in ROOT.rglob('*') if p.is_file() and (p.suffix in ['.py','.bat'] or p.name in ['RUN_CONFIG_V13_1.json','requirements.txt']) and not any(x in p.parts for x in ['.venv','__pycache__'])]
    code={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(active)}
    for p in sorted((ROOT/'raw_inputs').iterdir()): code[p.relative_to(ROOT).as_posix()]=sha(p)
    (ROOT/'FROZEN_CODE_HASHES_V13_1.json').write_text(json.dumps({'framework_version':'13.2.0-98','sha256':code},indent=2),encoding='utf-8')
    print('Frozen data, raw workbooks, active code and configuration. Formal results have not been generated.')
if __name__=='__main__': main()
