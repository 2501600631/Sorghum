from pathlib import Path
import json, tempfile
import pandas as pd
import numpy as np
from sorghum_v13_data_utils import load_data, _load_presentation_csv, make_regression_stratified_splits

ROOT=Path(__file__).resolve().parents[1]
def main():
    d=load_data(ROOT/'experimental_data',2)
    assert len(d.sample_ids)==98 and d.whole_replicates.shape==(98,1,1050)
    chemistry=pd.read_excel(ROOT/'raw_inputs/reference_chemistry_raw.xlsx').set_index('sample_id')
    np.testing.assert_allclose(d.chemistry.iloc[:,1:5].to_numpy(float),chemistry.loc[d.sample_ids].to_numpy(float),rtol=0,atol=1e-12)
    for sheet,actual in [('whole_grain_400_2499p5',d.whole_mean),('ground_grain_400_2499p5',d.ground_mean)]:
        raw=pd.read_excel(ROOT/'raw_inputs/paired_spectra_raw.xlsx',sheet_name=sheet).set_index('sample_id')
        np.testing.assert_allclose(actual,raw.loc[d.sample_ids].iloc[:,::4].to_numpy(float),rtol=0,atol=1e-12)

    assert '202512-06141' in d.sample_ids and '202512-06142' in d.sample_ids
    assert abs(d.chemistry.set_index('sample_id').loc['202512-06142','crude_protein']-9.20)<1e-12
    y=d.chemistry.crude_protein.to_numpy(float)
    splits=make_regression_stratified_splits(y,5,10,20260830)
    for repeat in range(10):
        tests=[]
        for r,f,tr,te in splits:
            if r!=repeat: continue
            assert not set(tr)&set(te) and len(tr)+len(te)==98
            tests.extend(te.tolist())
        assert sorted(tests)==list(range(98))

    groups=np.array([str(i//2) for i in range(98)])
    for _,_,tr,te in make_regression_stratified_splits(y,5,2,20260830,groups):
        assert not set(groups[tr])&set(groups[te])
    raw=pd.read_csv(ROOT/'experimental_data/whole_grain_replicate_spectra.csv.gz')
    rejects=[]
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp)/'bad.csv.gz'
        for case in ['duplicate_id','replicate_two','fractional_replicate','nonfinite']:
            bad=raw.copy()
            if case=='duplicate_id': bad.loc[1,'sample_id']=bad.loc[0,'sample_id']
            elif case=='replicate_two': bad.loc[0,'replicate']=2
            elif case=='fractional_replicate': bad.loc[0,'replicate']=1.5
            else: bad.iloc[0,2]=np.inf
            bad.to_csv(path,index=False,compression='gzip')
            try: _load_presentation_csv(path)
            except ValueError: rejects.append(case)
            else: raise AssertionError('Invalid input accepted: '+case)
    result={'status':'passed','samples':98,'source_values_match':True,'full_grid_decimation_matches':True,'outer_5x10_coverage_and_isolation':True,'group_isolation':True,'invalid_inputs_rejected':rejects}
    (ROOT/'docs/INDEPENDENT_SAMPLE_TEST.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
