from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
from typing import Dict, Iterable, Tuple
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import spearmanr

FRAMEWORK_VERSION="13.2.0-98"
DEFAULT_SEED=20260830
PUBLICATION_RANDOM_REPEATS=500
PUBLICATION_BOOTSTRAPS=2000
WHOLE_FAMILIES=["whole_pls","whole_ridge","whole_elastic","whole_ipls","whole_svr","whole_krr","whole_pls_krr","whole_extratrees"]
DEFAULT_FRACTIONS=[0.0,0.05,0.10,0.20,0.30,0.50,1.0]


def stable_seed(seed:int,*parts:object)->int:
    text="|".join([str(seed),*map(str,parts)])
    return int(hashlib.sha256(text.encode()).hexdigest()[:8],16)


def metrics(y:np.ndarray,p:np.ndarray)->Dict[str,float]:
    rmse=float(np.sqrt(mean_squared_error(y,p)))
    return {"r2":float(r2_score(y,p)),"rmse":rmse,"mae":float(mean_absolute_error(y,p)),"bias":float(np.mean(p-y)),"rpd":float(np.std(y,ddof=1)/rmse) if rmse>0 else np.inf}


def _cluster_col(df:pd.DataFrame)->str:
    return "inference_group" if "inference_group" in df.columns else "sample_id"


def paired_mse_delta(frame:pd.DataFrame,method_a:str,method_b:str)->pd.DataFrame:
    keys=["sample_id","outer_repeat","outer_fold"]
    if "inference_group" in frame.columns: keys.insert(1,"inference_group")
    a=frame[frame.method==method_a][[*keys,"target","y_true","y_pred"]].rename(columns={"y_pred":"pred_a"})
    b=frame[frame.method==method_b][[*keys,"target","y_pred"]].rename(columns={"y_pred":"pred_b"})
    m=a.merge(b,on=[*keys,"target"],validate="one_to_one")
    m["sqerr_a"]=(m.y_true-m.pred_a)**2; m["sqerr_b"]=(m.y_true-m.pred_b)**2
    m["mse_improvement_b_over_a"]=m.sqerr_a-m.sqerr_b
    return m


def cluster_bootstrap_delta(m:pd.DataFrame,seed:int,n_boot:int)->Tuple[float,float,float,str,int]:
    col=_cluster_col(m)
    per=m.groupby(col,sort=False)["mse_improvement_b_over_a"].mean().to_numpy(float)
    if len(per)<2: raise RuntimeError("At least two inference clusters are required")
    observed=float(per.mean()); rng=np.random.default_rng(seed); boots=np.empty(n_boot,float)
    for i in range(n_boot): boots[i]=float(rng.choice(per,size=len(per),replace=True).mean())
    lo,hi=np.quantile(boots,[.025,.975])
    return observed,float(lo),float(hi),col,int(len(per))


def utility_summary(oof:pd.DataFrame,seed:int,n_boot:int)->pd.DataFrame:
    rows=[]; targets=list(dict.fromkeys(oof.target.tolist()))
    comparisons=[("weak_baseline","whole_ridge","nested_ground_assisted"),("strong_baseline","nested_strong_whole","nested_ground_assisted"),("whole_vs_ground_reference","nested_strong_whole","ground_reference")]
    for t in targets:
        ft=oof[oof.target==t]
        for label,a,b in comparisons:
            if not {a,b}.issubset(set(ft.method)): continue
            m=paired_mse_delta(ft,a,b); obs,lo,hi,col,nc=cluster_bootstrap_delta(m,stable_seed(seed,t,label),n_boot)
            per_rep=m.groupby("outer_repeat").agg(mse_a=("sqerr_a","mean"),mse_b=("sqerr_b","mean"),delta=("mse_improvement_b_over_a","mean")).reset_index()
            rows.append({"target":t,"comparison":label,"method_a":a,"method_b":b,"positive_delta_means_better_b":True,
                         "observed_mean_mse_improvement_b_over_a":obs,"bootstrap95_low":lo,"bootstrap95_high":hi,
                         "inference_cluster_column":col,"n_inference_clusters":nc,
                         "median_repeat_mse_a":float(per_rep.mse_a.median()),"median_repeat_mse_b":float(per_rep.mse_b.median()),"median_repeat_delta":float(per_rep.delta.median()),
                         "interpretation":("evidence_of_incremental_value" if lo>0 else ("evidence_of_harm" if hi<0 else "no_clear_incremental_value"))})
    return pd.DataFrame(rows)


def build_policy_frame(ft:pd.DataFrame)->pd.DataFrame:
    keys=["sample_id","outer_repeat","outer_fold"]
    if "inference_group" in ft.columns: keys.insert(1,"inference_group")
    truth=ft[ft.method=="nested_strong_whole"][[*keys,"y_true","y_pred"]].rename(columns={"y_pred":"pred_whole"})
    ground=ft[ft.method=="ground_reference"][[*keys,"y_pred"]].rename(columns={"y_pred":"pred_ground"})
    wide=ft[ft.method.isin(WHOLE_FAMILIES)].pivot_table(index=keys,columns="method",values="y_pred",aggfunc="first")
    available=[c for c in WHOLE_FAMILIES if c in wide.columns]
    if len(available)<3: raise RuntimeError("Selective grinding requires at least three whole-family OOF streams")
    wide=wide.reset_index(); wide["uncertainty_score"]=wide[available].std(axis=1,ddof=1)
    out=truth.merge(ground,on=keys,validate="one_to_one").merge(wide[[*keys,"uncertainty_score"]],on=keys,validate="one_to_one")
    out["abs_whole_error"]=(out.y_true-out.pred_whole).abs()
    out["sqerr_benefit_ground_over_whole"]=(out.y_true-out.pred_whole)**2-(out.y_true-out.pred_ground)**2
    return out


def _allocate_fold_quotas(g:pd.DataFrame,fraction:float)->Dict[int,int]:
    sizes=g.groupby("outer_fold").size().sort_index()
    if fraction<=0: return {int(f):0 for f in sizes.index}
    if fraction>=1: return {int(f):int(n) for f,n in sizes.items()}
    target=int(round(len(g)*float(fraction))); raw=sizes.astype(float)*float(fraction); base=np.floor(raw).astype(int)
    remaining=target-int(base.sum()); rema=(raw-base).sort_values(ascending=False,kind="mergesort")
    quotas={int(f):int(base.loc[f]) for f in sizes.index}
    for f in list(rema.index)[:remaining]: quotas[int(f)]+=1
    if sum(quotas.values())!=target: raise RuntimeError("Fold quota allocation failed")
    return quotas


def _guided_indices(g:pd.DataFrame,quotas:Dict[int,int])->set:
    selected=set()
    for fold,gf in g.groupby("outer_fold",sort=True):
        q=int(quotas[int(fold)])
        if q: selected.update(gf.nlargest(q,"uncertainty_score").index.tolist())
    return selected


def evaluate_selective_grinding(oof:pd.DataFrame,seed:int,random_repeats:int,fractions:Iterable[float])->Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    policy_rows=[]; corr_rows=[]; sample_rows=[]
    for target in list(dict.fromkeys(oof.target.tolist())):
        ft=oof[oof.target==target]
        if not {"nested_strong_whole","ground_reference"}.issubset(set(ft.method)): continue
        pf=build_policy_frame(ft)
        for rep,g0 in pf.groupby("outer_repeat",sort=True):
            g=g0.copy().reset_index(drop=True)

            for fold,gf in g.groupby("outer_fold",sort=True):
                rho,pv=spearmanr(gf.uncertainty_score,gf.abs_whole_error)
                corr_rows.append({"target":target,"outer_repeat":int(rep),"outer_fold":int(fold),"n_test_samples":len(gf),"spearman_uncertainty_vs_abs_error":float(rho),"p_value_descriptive":float(pv)})
            n=len(g)
            for frac in fractions:
                quotas=_allocate_fold_quotas(g,float(frac)); k=sum(quotas.values()); idx_sel=_guided_indices(g,quotas)
                pred=g.pred_whole.to_numpy().copy()
                if idx_sel: pred[list(idx_sel)]=g.loc[list(idx_sel),"pred_ground"].to_numpy()
                met=metrics(g.y_true.to_numpy(),pred)
                policy_rows.append({"target":target,"outer_repeat":int(rep),"policy":"uncertainty_guided","requested_fraction":float(frac),"n_ground":k,"grinding_burden":k/n,"selection_scope":"within_outer_fold",**met})
                for i,row in g.iterrows():
                    sample_rows.append({"target":target,"outer_repeat":int(rep),"outer_fold":int(row.outer_fold),"sample_id":row.sample_id,"inference_group":getattr(row,"inference_group",row.sample_id),"policy":"uncertainty_guided","requested_fraction":float(frac),"selected_for_grinding":i in idx_sel,"fold_quota":quotas[int(row.outer_fold)],"uncertainty_score":row.uncertainty_score,"y_true":row.y_true,"y_pred":pred[i]})
                random_metrics=[]
                for rr in range(random_repeats):
                    rng=np.random.default_rng(stable_seed(seed,target,rep,frac,rr,"random")); ridx=set()
                    for fold,gf in g.groupby("outer_fold",sort=True):
                        q=quotas[int(fold)]
                        if q: ridx.update(rng.choice(gf.index.to_numpy(int),size=q,replace=False).tolist())
                    rp=g.pred_whole.to_numpy().copy()
                    if ridx: rp[list(ridx)]=g.loc[list(ridx),"pred_ground"].to_numpy()
                    random_metrics.append(metrics(g.y_true.to_numpy(),rp))
                rmean={m:float(np.mean([x[m] for x in random_metrics])) for m in ["r2","rmse","mae","bias","rpd"]}
                rq={m:(float(np.quantile([x[m] for x in random_metrics],.025)),float(np.quantile([x[m] for x in random_metrics],.975))) for m in ["r2","rmse","mae","bias","rpd"]}
                policy_rows.append({"target":target,"outer_repeat":int(rep),"policy":"random_matched","requested_fraction":float(frac),"n_ground":k,"grinding_burden":k/n,"selection_scope":"within_outer_fold",**rmean,"random_r2_q025":rq["r2"][0],"random_r2_q975":rq["r2"][1],"random_rmse_q025":rq["rmse"][0],"random_rmse_q975":rq["rmse"][1]})

                oid=set()
                for fold,gf in g.groupby("outer_fold",sort=True):
                    q=quotas[int(fold)]
                    if q: oid.update(gf.nlargest(q,"sqerr_benefit_ground_over_whole").index.tolist())
                op=g.pred_whole.to_numpy().copy()
                if oid: op[list(oid)]=g.loc[list(oid),"pred_ground"].to_numpy()
                policy_rows.append({"target":target,"outer_repeat":int(rep),"policy":"oracle_upper_bound","requested_fraction":float(frac),"n_ground":k,"grinding_burden":k/n,"selection_scope":"within_outer_fold",**metrics(g.y_true.to_numpy(),op)})
    return pd.DataFrame(policy_rows),pd.DataFrame(corr_rows),pd.DataFrame(sample_rows)


def summarize_policies(pol:pd.DataFrame)->pd.DataFrame:
    rows=[]
    for keys,g in pol.groupby(["target","policy","requested_fraction"],sort=False):
        row={"target":keys[0],"policy":keys[1],"requested_fraction":keys[2],"median_grinding_burden":float(g.grinding_burden.median())}
        for m in ["r2","rmse","mae","bias","rpd"]:
            row[f"median_{m}"]=float(g[m].median()); row[f"q25_{m}"]=float(g[m].quantile(.25)); row[f"q75_{m}"]=float(g[m].quantile(.75))
        rows.append(row)
    return pd.DataFrame(rows)


def enforce_lock(args:argparse.Namespace)->None:
    if args.profile!="publication": return
    expected={"seed":DEFAULT_SEED,"random_repeats":PUBLICATION_RANDOM_REPEATS,"bootstrap_repeats":PUBLICATION_BOOTSTRAPS}
    bad=[f"{k}={getattr(args,k)!r} (locked={v!r})" for k,v in expected.items() if getattr(args,k)!=v]
    if [float(x) for x in args.fractions]!=DEFAULT_FRACTIONS: bad.append(f"fractions={args.fractions!r} (locked={DEFAULT_FRACTIONS!r})")
    if bad: raise RuntimeError("V13.2.0-98 UTILITY DESIGN LOCK VIOLATION: "+"; ".join(bad))


def build_parser()->argparse.ArgumentParser:
    ap=argparse.ArgumentParser(); root=Path(__file__).resolve().parents[1]
    ap.add_argument("--results-dir",type=Path,default=root/"v13_2_results"); ap.add_argument("--profile",choices=["publication","smoke"],default="publication")
    ap.add_argument("--seed",type=int,default=DEFAULT_SEED); ap.add_argument("--random-repeats",type=int,default=PUBLICATION_RANDOM_REPEATS); ap.add_argument("--bootstrap-repeats",type=int,default=PUBLICATION_BOOTSTRAPS); ap.add_argument("--fractions",nargs="+",type=float,default=DEFAULT_FRACTIONS)
    return ap

def main()->None:
    args=build_parser().parse_args(); enforce_lock(args)
    pred=args.results_dir/"01_predictions"/"all_targets_outer_oof_predictions.csv"
    if not pred.exists(): raise FileNotFoundError(pred)
    oof=pd.read_csv(pred); out=args.results_dir/"05_utility_boundary"; out.mkdir(parents=True,exist_ok=True)
    utility_summary(oof,args.seed,args.bootstrap_repeats).to_csv(out/"trait_privileged_incremental_value.csv",index=False)
    policies,corr,sample=evaluate_selective_grinding(oof,args.seed,args.random_repeats,args.fractions)
    policies.to_csv(out/"selective_grinding_repeat_performance.csv",index=False); summarize_policies(policies).to_csv(out/"selective_grinding_cost_performance_summary.csv",index=False); corr.to_csv(out/"uncertainty_validity_summary.csv",index=False); sample.to_csv(out/"selective_grinding_sample_decisions.csv.gz",index=False,compression="gzip")
    summary={"framework_version":FRAMEWORK_VERSION,"status":"completed","targets":list(dict.fromkeys(oof.target.tolist())),"selective_grinding_is_exploratory":True,"uncertainty_policy_uses_y_true_for_selection":False,"uncertainty_definition":"SD across fixed whole-only family OOF predictions, ranked only among test samples produced by the same outer-fold fitted model","selection_scope":"within_outer_fold","oracle_uses_y_true":True,"oracle_is_deployable":False,"random_policy_repeats":args.random_repeats,"fractions":args.fractions,"important_interpretation":"Retrospective OOF decision-policy simulation; not independent deployment validation."}
    (out/"UTILITY_BOUNDARY_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print("V13.2.0-98 utility-boundary analysis complete:",out)

if __name__=="__main__": main()
