from pathlib import Path
import json
import numpy as np
import sorghum_v13_1_frozen as c
import v13_1_privileged_utility_analysis as u


def main():
    root=Path(__file__).resolve().parents[1]
    D=c.load_data(root/'experimental_data',2.0,400.0,2498.0)
    y=D.chemistry.set_index('sample_id').loc[D.sample_ids,'crude_protein'].to_numpy(float)
    n=len(D.sample_ids); cut=int(n*0.7); tr=np.arange(cut); te=np.arange(cut,n)
    reps=c.representation_space('publication'); rep=next(r for r in reps if r.name.startswith('FULL__SNV'))
    results=[]
    groups=[('whole',c._candidate_grids(rep,'publication',False)),('ground',c._ground_candidate_grids(rep,'publication')),('transfer',c._candidate_grids(rep,'publication',True))]
    for group,cands in groups:
        fams={}
        for cand in cands: fams.setdefault(cand.family,cand)
        for fam,cand in sorted(fams.items()):
            if group=='ground':
                fit=c.fit_ground_candidate(cand,rep,D.wavelengths,D.ground_mean[tr],y[tr],12345)
                pred=c.predict_ground_candidate(fit,D.ground_mean[te])
            else:
                fit=c.fit_candidate(cand,rep,D.wavelengths,D.whole_mean[tr],D.ground_mean[tr],y[tr],12345)
                pred=c.predict_candidate(fit,D.whole_mean[te])
            ok=bool(len(pred)==len(te) and np.isfinite(pred).all())
            results.append({'group':group,'family':fam,'candidate_key':cand.key,'finite_predictions':ok,'n_test':len(pred)})
            if not ok: raise RuntimeError(f'Family self-test failed: {group}/{fam}')


    rng=np.random.default_rng(20260902)
    Xw_syn=rng.normal(size=(18,12))
    shift=np.linspace(0.2,1.4,12)
    Xg_syn=Xw_syn+shift
    y_syn=np.linspace(0.0,1.0,18)
    epo=c.fit_epo_base(Xw_syn,Xg_syn,y_syn,1,'ridge',{'alpha':1.0},123)
    v=epo['V'][:,0]; shift_u=shift/np.linalg.norm(shift)
    epo_alignment=float(abs(np.dot(v,shift_u)))
    if epo_alignment < 0.999999:
        raise RuntimeError(f'Non-centered EPO self-test failed: alignment={epo_alignment}')



    import pandas as pd
    probe=pd.DataFrame([
        {'candidate_key':'lower_rmse_high_complexity','mean_rmse':0.500,'se_rmse':0.100,'complexity':5},
        {'candidate_key':'higher_rmse_low_complexity','mean_rmse':0.520,'se_rmse':0.001,'complexity':1},
    ])
    selected=str(c.select_min_inner_rmse(probe)['candidate_key'])
    if selected!='lower_rmse_high_complexity':
        raise RuntimeError(f'Minimum-inner-RMSE selector self-test failed: {selected}')



    import pandas as pd
    quota_probe=pd.DataFrame({'outer_fold':[i % 5 for i in range(n)]})
    q10=u._allocate_fold_quotas(quota_probe,0.10); q20=u._allocate_fold_quotas(quota_probe,0.20)
    if sum(q10.values())!=round(n*0.10) or sum(q20.values())!=round(n*0.20):
        raise RuntimeError(f'Fold-local quota self-test failed: {q10}, {q20}')

    report={'status':'passed','framework_version':'13.2.0-98','purpose':'Executable QA only; fixed 70-percent/train remaining/test split is not scientific validation.',
            'fold_local_quota_10pct':q10,'fold_local_quota_20pct':q20,
            'representation':rep.name,'families_tested':len(results),'epo_noncentered_shift_alignment':epo_alignment,
            'candidate_selection_selftest':selected,'results':results}
    out=root/'docs'/'V13_1_ALGORITHM_FAMILY_SELFTEST.json'; out.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'status':'passed','families_tested':len(results),'output':str(out)},indent=2))
if __name__=='__main__': main()
