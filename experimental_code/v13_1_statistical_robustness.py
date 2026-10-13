from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
from typing import Dict, Iterable, Sequence, Tuple
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from sklearn.metrics import r2_score
from sorghum_v13_data_utils import load_data

FRAMEWORK_VERSION="13.2.0-98"
DEFAULT_SEED=20260830
PRIMARY_TARGET="crude_protein"
WHOLE_FAMILIES=["whole_pls","whole_ridge","whole_elastic","whole_ipls","whole_svr","whole_krr","whole_pls_krr","whole_extratrees"]
PUBLICATION_MARGINS_PCT=(2.5,5.0,10.0)
PUBLICATION_BUDGETS=(0.10,0.20)
PUBLICATION_BOOTSTRAPS=2000
PUBLICATION_SIGNFLIPS=20000
PUBLICATION_LABEL_PERMUTATIONS=20000


def stable_seed(seed:int,*parts:object)->int:
    return int(hashlib.sha256("|".join([str(seed),*map(str,parts)]).encode()).hexdigest()[:8],16)


def holm(pvals:Sequence[float])->np.ndarray:
    p=np.asarray(pvals,float); order=np.argsort(p); out=np.empty_like(p); running=0.0; m=len(p)
    for rank,idx in enumerate(order):
        val=min(1.0,(m-rank)*p[idx]); running=max(running,val); out[idx]=running
    return out


def _cluster_col(df:pd.DataFrame)->str:
    return "inference_group" if "inference_group" in df.columns else "sample_id"


def _pair_methods(oof:pd.DataFrame,target:str,a:str,b:str)->pd.DataFrame:
    keys=["sample_id","outer_repeat","outer_fold"]
    if "inference_group" in oof.columns: keys.insert(1,"inference_group")
    aa=oof[(oof.target==target)&(oof.method==a)][[*keys,"y_true","y_pred"]].rename(columns={"y_pred":"pred_a"})
    bb=oof[(oof.target==target)&(oof.method==b)][[*keys,"y_pred"]].rename(columns={"y_pred":"pred_b"})
    m=aa.merge(bb,on=keys,validate="one_to_one")
    if m.empty: raise RuntimeError(f"No paired OOF rows for {target}: {a} vs {b}")
    return m


def _rmse(y,p)->float:
    y=np.asarray(y,float); p=np.asarray(p,float); return float(np.sqrt(np.mean((y-p)**2)))


def _cluster_boot_rmse_gain(m:pd.DataFrame,seed:int,n_boot:int)->Dict[str,float]:
    col=_cluster_col(m); clusters=np.array(sorted(m[col].astype(str).unique()),dtype=str)
    rmse_a=_rmse(m.y_true,m.pred_a); rmse_b=_rmse(m.y_true,m.pred_b); delta=rmse_a-rmse_b; rel=100.0*delta/rmse_a if rmse_a>0 else np.nan
    by={c:m[m[col].astype(str)==c] for c in clusters}; rng=np.random.default_rng(seed); bd=np.empty(n_boot); br=np.empty(n_boot)
    for i in range(n_boot):
        chosen=rng.choice(clusters,size=len(clusters),replace=True); ys=[]; pa=[]; pb=[]
        for c in chosen:
            g=by[str(c)]; ys.append(g.y_true.to_numpy(float)); pa.append(g.pred_a.to_numpy(float)); pb.append(g.pred_b.to_numpy(float))
        y=np.concatenate(ys); a=np.concatenate(pa); b=np.concatenate(pb); ra=_rmse(y,a); rb=_rmse(y,b); bd[i]=ra-rb; br[i]=100.0*(ra-rb)/ra if ra>0 else np.nan
    qd=np.quantile(bd,[.025,.975]); qr=np.nanquantile(br,[.025,.975])
    return {"rmse_a":rmse_a,"rmse_b":rmse_b,"delta_rmse_b_better_positive":delta,"relative_rmse_improvement_pct_b_over_a":rel,
            "delta_rmse_bootstrap95_low":float(qd[0]),"delta_rmse_bootstrap95_high":float(qd[1]),"relative_improvement_bootstrap95_low_pct":float(qr[0]),"relative_improvement_bootstrap95_high_pct":float(qr[1]),
            "inference_cluster_column":col,"n_inference_clusters":int(len(clusters))}


def practical_gain_sensitivity(oof:pd.DataFrame,seed:int,n_boot:int,margins:Iterable[float])->pd.DataFrame:
    rows=[]
    for target in list(dict.fromkeys(oof.target.tolist())):
        if not {"nested_strong_whole","nested_ground_assisted"}.issubset(set(oof[oof.target==target].method)): continue
        m=_pair_methods(oof,target,"nested_strong_whole","nested_ground_assisted"); base=_cluster_boot_rmse_gain(m,stable_seed(seed,"practical_gain",target),n_boot)
        for margin in margins:
            lo=base["relative_improvement_bootstrap95_low_pct"]; hi=base["relative_improvement_bootstrap95_high_pct"]
            interp="gain_of_at_least_margin_excluded_by_95pct_CI" if hi<float(margin) else ("gain_exceeding_margin_supported_by_95pct_CI" if lo>float(margin) else "margin_crossed_by_95pct_CI")
            rows.append({"target":target,"method_a":"nested_strong_whole","method_b":"nested_ground_assisted",**base,"sensitivity_margin_pct":float(margin),"margin_is_reference_method_or_industry_derived":False,"interpretation":interp,
                         "guardrail":"Prespecified sensitivity threshold, not an assay-derived equivalence margin."})
    return pd.DataFrame(rows)


def build_policy_frame(oof:pd.DataFrame,target:str)->pd.DataFrame:
    ft=oof[oof.target==target]; keys=["sample_id","outer_repeat","outer_fold"]
    if "inference_group" in ft.columns: keys.insert(1,"inference_group")
    truth=ft[ft.method=="nested_strong_whole"][[*keys,"y_true","y_pred"]].rename(columns={"y_pred":"pred_whole"})
    ground=ft[ft.method=="ground_reference"][[*keys,"y_pred"]].rename(columns={"y_pred":"pred_ground"})
    wide=ft[ft.method.isin(WHOLE_FAMILIES)].pivot_table(index=keys,columns="method",values="y_pred",aggfunc="first")
    available=[c for c in WHOLE_FAMILIES if c in wide.columns]
    if len(available)<3: raise RuntimeError("At least three whole-family OOF streams are required")
    wide=wide.reset_index(); wide["uncertainty_score"]=wide[available].std(axis=1,ddof=1)
    m=truth.merge(ground,on=keys,validate="one_to_one").merge(wide[[*keys,"uncertainty_score"]],on=keys,validate="one_to_one")
    m["sqerr_whole"]=(m.y_true-m.pred_whole)**2; m["sqerr_ground"]=(m.y_true-m.pred_ground)**2
    return m


def _allocate_fold_quotas(g:pd.DataFrame,budget:float)->Dict[int,int]:
    sizes=g.groupby("outer_fold").size().sort_index(); target=int(round(len(g)*float(budget))); raw=sizes.astype(float)*float(budget); base=np.floor(raw).astype(int)
    rem=target-int(base.sum()); frac=(raw-base).sort_values(ascending=False,kind="mergesort"); q={int(f):int(base.loc[f]) for f in sizes.index}
    for f in list(frac.index)[:rem]: q[int(f)]+=1
    return q


def _signflip_cluster_means(vals:np.ndarray,seed:int,n_perm:int)->Tuple[float,float]:
    vals=np.asarray(vals,float); obs=float(vals.mean()); rng=np.random.default_rng(seed); null=np.empty(n_perm)
    for i in range(n_perm): null[i]=float(np.mean(vals*rng.choice([-1.0,1.0],size=len(vals))))
    return obs,float((1+np.sum(null>=obs))/(n_perm+1))


def selective_grinding_inference(oof:pd.DataFrame,seed:int,n_boot:int,n_perm:int,budgets:Iterable[float])->pd.DataFrame:
    rows=[]
    for target in list(dict.fromkeys(oof.target.tolist())):
        pf=build_policy_frame(oof,target); col=_cluster_col(pf)
        for budget in budgets:
            parts=[]
            for rep,g0 in pf.groupby("outer_repeat",sort=True):
                g=g0.copy().reset_index(drop=True); quotas=_allocate_fold_quotas(g,float(budget)); selected=set()
                fold_fraction={}
                for fold,gf in g.groupby("outer_fold",sort=True):
                    q=quotas[int(fold)]; fold_fraction[int(fold)]=q/len(gf)
                    if q: selected.update(gf.nlargest(q,"uncertainty_score").index.tolist())
                g["selected_guided"]=[i in selected for i in range(len(g))]
                g["sqerr_guided"]=np.where(g.selected_guided,g.sqerr_ground,g.sqerr_whole)
                g["random_selection_probability"]=[fold_fraction[int(f)] for f in g.outer_fold]
                g["sqerr_random_expected"]=(1.0-g.random_selection_probability)*g.sqerr_whole+g.random_selection_probability*g.sqerr_ground
                g["delta_sqerr_guided_better_positive"]=g.sqerr_random_expected-g.sqerr_guided
                g["actual_fraction"]=sum(quotas.values())/len(g); parts.append(g)
            m=pd.concat(parts,ignore_index=True); per_cluster=m.groupby(col,sort=False)["delta_sqerr_guided_better_positive"].mean().to_numpy(float)
            obs_mse,p=_signflip_cluster_means(per_cluster,stable_seed(seed,"selective_signflip",target,budget),n_perm)
            clusters=np.array(sorted(m[col].astype(str).unique()),dtype=str); by={c:m[m[col].astype(str)==c] for c in clusters}; rng=np.random.default_rng(stable_seed(seed,"selective_boot",target,budget)); bdm=np.empty(n_boot); bdr=np.empty(n_boot)
            for i in range(n_boot):
                chosen=rng.choice(clusters,size=len(clusters),replace=True); rg=[]; rr=[]
                for c in chosen:
                    gg=by[str(c)]; rg.append(gg.sqerr_guided.to_numpy(float)); rr.append(gg.sqerr_random_expected.to_numpy(float))
                eg=np.concatenate(rg); er=np.concatenate(rr); bdm[i]=float(np.mean(er-eg)); bdr[i]=float(np.sqrt(np.mean(er))-np.sqrt(np.mean(eg)))
            ci_m=np.quantile(bdm,[.025,.975]); ci_r=np.quantile(bdr,[.025,.975]); mse_g=float(m.sqerr_guided.mean()); mse_r=float(m.sqerr_random_expected.mean())
            rows.append({"target":target,"requested_budget":float(budget),"median_actual_fraction":float(m.groupby("outer_repeat").actual_fraction.first().median()),"selection_scope":"within_outer_fold",
                         "guided_rmse":math.sqrt(mse_g),"matched_random_expected_rmse":math.sqrt(mse_r),"delta_rmse_guided_better_positive":math.sqrt(mse_r)-math.sqrt(mse_g),"delta_mse_guided_better_positive":obs_mse,
                         "delta_mse_bootstrap95_low":float(ci_m[0]),"delta_mse_bootstrap95_high":float(ci_m[1]),"delta_rmse_bootstrap95_low":float(ci_r[0]),"delta_rmse_bootstrap95_high":float(ci_r[1]),
                         "one_sided_cluster_signflip_p":p,"inference_cluster_column":col,"n_inference_clusters":len(clusters),"signflip_permutations":int(n_perm),"bootstrap_repeats":int(n_boot),
                         "random_comparator_definition":"Exact expected squared error under uniform selection of the prespecified quota within each outer fold.",
                         "bootstrap_conditioning":"Conditional on frozen OOF uncertainty scores and fold-local guided decisions; inference_group is the resampling cluster when available."})
    df=pd.DataFrame(rows)
    if len(df): df["holm_adjusted_p_all_target_budget_tests"]=holm(df["one_sided_cluster_signflip_p"].to_numpy(float))
    return df


def corrected_repeated_cv_sensitivity(oof:pd.DataFrame)->pd.DataFrame:






    rows=[]
    for target in list(dict.fromkeys(oof.target.tolist())):
        m=_pair_methods(oof,target,"nested_strong_whole","nested_ground_assisted"); m["delta_sqerr"]=(m.y_true-m.pred_a)**2-(m.y_true-m.pred_b)**2
        fold=m.groupby(["outer_repeat","outer_fold"]).agg(delta_mse=("delta_sqerr","mean"),n_test=("sample_id","size")).reset_index(); total=int(m.sample_id.nunique()); fold["n_train"]=total-fold.n_test
        d=fold.delta_mse.to_numpy(float); n=len(d); mean=float(d.mean()); var=float(np.var(d,ddof=1)); ratio=float(np.mean(fold.n_test/fold.n_train)); se=float(np.sqrt((1.0/n+ratio)*var)) if var>0 else 0.0
        tstat=float(mean/se) if se>0 else (float("inf") if mean>0 else float("-inf") if mean<0 else 0.0); p=float(student_t.sf(tstat,df=n-1)) if np.isfinite(tstat) else (0.0 if tstat>0 else 1.0)
        rows.append({"target":target,"method_a":"nested_strong_whole","method_b":"nested_ground_assisted","positive_delta_means_better_b":True,"n_outer_fold_differences":n,"mean_fold_mse_improvement_b_over_a":mean,
                     "fold_difference_sd":float(np.sqrt(var)),"mean_test_train_ratio":ratio,"corrected_standard_error":se,"corrected_t":tstat,"df":n-1,"one_sided_p":p,
                     "test_definition":"Nadeau-Bengio-style corrected repeated-CV t sensitivity analysis on paired outer-fold MSE differences (5x10 in publication mode).","interpretation_guardrail":"Supplemental development-cohort sensitivity analysis; not independent confirmation."})
    return pd.DataFrame(rows)


def primary_oof_label_permutation(oof:pd.DataFrame,seed:int,n_perm:int)->Tuple[pd.DataFrame,Dict[str,object]]:
    g=oof[(oof.target==PRIMARY_TARGET)&(oof.method=="nested_strong_whole")].copy(); col=_cluster_col(g)
    if col!="sample_id" and g.groupby(col).sample_id.nunique().max()>1:
        return pd.DataFrame(columns=["permutation","null_r2","null_rmse"]),{"status":"skipped_grouped_provenance","target":PRIMARY_TARGET,"method":"nested_strong_whole","inference_cluster_column":col,"reason":"Sample-level outcome permutation is not valid when multiple samples share an inference group; test intentionally omitted.","not_full_y_randomization":True,"permutations":0}
    sample_y=g.groupby("sample_id").y_true.first().sort_index(); sids=sample_y.index.to_numpy(str); yvals=sample_y.to_numpy(float); pred_by_sid={sid:g[g.sample_id==sid].sort_values("outer_repeat").y_pred.to_numpy(float) for sid in sids}
    pred=np.concatenate([pred_by_sid[sid] for sid in sids]); y=np.concatenate([np.repeat(yvals[i],len(pred_by_sid[sid])) for i,sid in enumerate(sids)]); obs_r2=float(r2_score(y,pred)); obs_rmse=_rmse(y,pred); rng=np.random.default_rng(seed); null_r2=np.empty(n_perm); null_rmse=np.empty(n_perm); reps=np.array([len(pred_by_sid[sid]) for sid in sids],int)
    for i in range(n_perm):
        perm=rng.permutation(yvals); yp=np.concatenate([np.repeat(perm[j],reps[j]) for j in range(len(sids))]); null_r2[i]=float(r2_score(yp,pred)); null_rmse[i]=_rmse(yp,pred)
    dist=pd.DataFrame({"permutation":np.arange(1,n_perm+1),"null_r2":null_r2,"null_rmse":null_rmse}); summary={"status":"completed","target":PRIMARY_TARGET,"method":"nested_strong_whole","observed_pooled_oof_r2":obs_r2,"observed_pooled_oof_rmse":obs_rmse,"permutations":int(n_perm),"sample_clusters":len(sids),"one_sided_p_r2_ge_observed":float((1+np.sum(null_r2>=obs_r2))/(n_perm+1)),"one_sided_p_rmse_le_observed":float((1+np.sum(null_rmse<=obs_rmse))/(n_perm+1)),"test_definition":"Biological-sample OOF prediction-label permutation; frozen predictions are not retrained.","not_full_y_randomization":True}
    return dist,summary


def build_parser()->argparse.ArgumentParser:
    ap=argparse.ArgumentParser(); root=Path(__file__).resolve().parents[1]; ap.add_argument("--results-dir",type=Path,default=root/"v13_2_results"); ap.add_argument("--data-dir",type=Path,default=root/"experimental_data"); ap.add_argument("--profile",choices=["publication","smoke"],default="publication"); ap.add_argument("--seed",type=int,default=DEFAULT_SEED); ap.add_argument("--bootstrap-repeats",type=int,default=PUBLICATION_BOOTSTRAPS); ap.add_argument("--signflip-repeats",type=int,default=PUBLICATION_SIGNFLIPS); ap.add_argument("--label-permutations",type=int,default=PUBLICATION_LABEL_PERMUTATIONS); return ap


def enforce_lock(args:argparse.Namespace)->None:
    if args.profile!="publication": return
    expected={"seed":DEFAULT_SEED,"bootstrap_repeats":PUBLICATION_BOOTSTRAPS,"signflip_repeats":PUBLICATION_SIGNFLIPS,"label_permutations":PUBLICATION_LABEL_PERMUTATIONS}; bad=[f"{k}={getattr(args,k)} (locked={v})" for k,v in expected.items() if getattr(args,k)!=v]
    if bad: raise RuntimeError("V13.2.0-98 ROBUSTNESS DESIGN LOCK VIOLATION: "+"; ".join(bad))


def main()->None:
    args=build_parser().parse_args(); enforce_lock(args); pred=args.results_dir/"01_predictions"/"all_targets_outer_oof_predictions.csv"
    if not pred.exists(): raise FileNotFoundError(pred)
    oof=pd.read_csv(pred); out=args.results_dir/"06_robustness"; out.mkdir(parents=True,exist_ok=True)
    practical_gain_sensitivity(oof,args.seed,args.bootstrap_repeats,PUBLICATION_MARGINS_PCT).to_csv(out/"practical_gain_sensitivity.csv",index=False)
    selective_grinding_inference(oof,args.seed,args.bootstrap_repeats,args.signflip_repeats,PUBLICATION_BUDGETS).to_csv(out/"selective_grinding_vs_random_inference.csv",index=False)
    corrected_repeated_cv_sensitivity(oof).to_csv(out/"corrected_repeated_cv_sensitivity.csv",index=False)
    (out/"technical_replicate_status.json").write_text(json.dumps({"status":"not_applicable","observations_per_sample_per_presentation":1,"reason":"98 independent samples; suffixes 1/2 are not technical replicates."},indent=2),encoding="utf-8")
    dist,perm_summary=primary_oof_label_permutation(oof,stable_seed(args.seed,"primary_label_permutation"),args.label_permutations); dist.to_csv(out/"primary_protein_oof_label_permutation_null.csv.gz",index=False,compression="gzip"); (out/"primary_protein_oof_label_permutation_summary.json").write_text(json.dumps(perm_summary,indent=2),encoding="utf-8")
    summary={"framework_version":FRAMEWORK_VERSION,"status":"completed","practical_gain_sensitivity_margins_pct":list(PUBLICATION_MARGINS_PCT),"selective_grinding_primary_exploratory_budgets":list(PUBLICATION_BUDGETS),"selective_grinding_selection_scope":"within_outer_fold","sample_or_source_group_cluster_bootstrap_repeats":args.bootstrap_repeats,"cluster_signflip_repeats":args.signflip_repeats,"primary_oof_label_permutations_requested":args.label_permutations,"corrected_repeated_cv_sensitivity":True,"external_validation_performed":False,"guardrails":["Practical-gain margins are sensitivity thresholds, not assay-derived equivalence margins.","Selective-grinding inference is retrospective OOF policy evidence, not deployment validation.","Corrected repeated-CV inference is supplemental development-cohort sensitivity evidence.","OOF label permutation does not retrain the full pipeline and is not full Y-randomization."]}
    (out/"ROBUSTNESS_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8"); print("V13.2.0-98 statistical robustness complete:",out)

if __name__=="__main__": main()
