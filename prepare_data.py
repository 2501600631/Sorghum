from pathlib import Path
import hashlib, json
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parent
TARGETS=['moisture','crude_protein','crude_ash','crude_fiber']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
    out=ROOT/'experimental_data'; out.mkdir(exist_ok=True)
    raw=ROOT/'raw_inputs'
    chem=pd.read_excel(raw/'reference_chemistry_raw.xlsx')
    if chem.shape!=(98,5): raise ValueError('Expected 98 chemistry rows and five columns')
    if list(chem.columns[1:])!=['moisture_pct','crude_protein_pct','crude_ash_pct','crude_fiber_pct']:
        raise ValueError('Unexpected chemistry column names/order')
    chem.columns=['sample_id']+TARGETS
    chem.sample_id=chem.sample_id.astype(str).str.strip()
    if chem.sample_id.duplicated().any() or chem.sample_id.eq('').any(): raise ValueError('Invalid chemistry IDs')
    if not np.isfinite(chem[TARGETS].to_numpy(float)).all(): raise ValueError('Invalid chemistry values')
    book=pd.ExcelFile(raw/'paired_spectra_raw.xlsx')
    sheets={'whole':'whole_grain_400_2499p5','ground':'ground_grain_400_2499p5'}
    audit={'independent_samples':98,'independence_basis':'User declaration on 2026-09-15',
           'id_rule':'Retain full IDs including final 1/2; no cross-sample averaging',
           'technical_replicates_available':False,'source_sha256':{p.name:sha(p) for p in raw.iterdir() if p.is_file()},'presentations':{}}
    mapping=chem[['sample_id']].copy(); mapping['chemistry_excel_row']=np.arange(2,100)
    for presentation,sheet in sheets.items():
        df=pd.read_excel(book,sheet_name=sheet)
        ids=df.iloc[:,0].astype(str).str.strip()
        if df.shape!=(98,4201) or ids.duplicated().any() or set(ids)!=set(chem.sample_id):
            raise ValueError('Spectrum count, duplicate IDs or sample alignment failure: '+sheet)
        waves=np.asarray(df.columns[1:],float)
        if not np.array_equal(waves,np.arange(400,2500,.5)): raise ValueError('Unexpected wavelength grid')
        vals=df.iloc[:,1:].to_numpy(float)
        if not np.isfinite(vals).all(): raise ValueError('Non-finite spectra')
        prepared=pd.DataFrame(vals,columns=['wl_'+str(w) for w in waves])
        prepared.insert(0,'replicate',1); prepared.insert(0,'sample_id',ids)
        prepared.to_csv(out/(presentation+'_grain_replicate_spectra.csv.gz'),index=False,compression={'method':'gzip','mtime':0})
        rowmap=dict(zip(ids,range(2,100))); mapping[presentation+'_excel_row']=mapping.sample_id.map(rowmap)
        audit['presentations'][presentation]={'source_sheet':sheet,'rows':98,'wavelengths':4200,'range_nm':[400,2499.5],'duplicate_spectral_rows':int(pd.DataFrame(vals).duplicated().sum())}
    chem.to_csv(out/'reference_chemistry.csv',index=False)
    mapping.to_csv(out/'sample_alignment.csv',index=False)

    mp=out/'sample_metadata_template.csv'
    if not mp.exists():
        meta=chem[['sample_id']].copy()
        meta['independence_status']='user_declared_independent_2026-09-15'
        for col in ['source_group','origin','batch','cultivar','laboratory_record']: meta[col]=''
        meta.to_csv(mp,index=False)
    rp=out/'reference_method_metadata_template.csv'
    if not rp.exists():
        ref=pd.DataFrame({'target':TARGETS})
        for col in ['standard_method_or_protocol','analytical_replicates','reporting_basis','reference_instrument_or_apparatus','laboratory_or_institution']: ref[col]=''
        ref.to_csv(rp,index=False)
    audit['chemistry_summary']=chem[TARGETS].describe().to_dict()
    (out/'DATA_IMPORT_AUDIT.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
    print('PASS: 98 unique independent IDs, matched whole/ground spectra and four chemistry targets; no averaging.')
if __name__=='__main__': prepare()
