#!/usr/bin/env python3
"""Data loading, grouped splitting, metrics, and reproducibility utilities for V13.

This module contains no model search or result-driven optimization. It defines the statistical
unit, paired spectral input checks, frozen-grid decimation, deterministic grouped regression
splits, and common output utilities used by the fairness-corrected V13 workflow.
"""
from __future__ import annotations
import hashlib, json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

TARGETS = ["moisture", "crude_protein", "crude_ash", "crude_fiber"]
TARGET_LABELS = {
    "moisture": "Moisture",
    "crude_protein": "Crude protein",
    "crude_ash": "Crude ash",
    "crude_fiber": "Crude fiber",
}


def stable_seed(base: int, *parts: Any) -> int:
    text = "|".join([str(base)] + [str(x) for x in parts])
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16) % 2_000_000_000


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b: break
            h.update(b)
    return h.hexdigest()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_json(obj: Any, path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def regression_metrics(y: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    y=np.asarray(y,float); pred=np.asarray(pred,float)
    rmse=float(np.sqrt(mean_squared_error(y,pred)))
    sd=float(np.std(y,ddof=1)) if len(y)>1 else float("nan")
    return {
        "r2":float(r2_score(y,pred)), "rmse":rmse,
        "mae":float(mean_absolute_error(y,pred)), "bias":float(np.mean(pred-y)),
        "rpd":float(sd/rmse) if rmse>0 and np.isfinite(sd) else float("nan"),
    }


def parse_wavelength_column(name: str) -> Optional[float]:
    if not name.startswith("wl_"): return None
    try: return float(name[3:])
    except ValueError: return None


@dataclass
class PairedSpectralData:
    sample_ids: np.ndarray
    wavelengths: np.ndarray
    whole_replicates: np.ndarray
    ground_replicates: np.ndarray
    chemistry: pd.DataFrame
    @property
    def whole_mean(self) -> np.ndarray: return self.whole_replicates.mean(axis=1)
    @property
    def ground_mean(self) -> np.ndarray: return self.ground_replicates.mean(axis=1)


def _load_presentation_csv(path: Path) -> Tuple[np.ndarray,np.ndarray,np.ndarray]:
    df=pd.read_csv(path,compression="infer")
    if not {"sample_id","replicate"}.issubset(df.columns):
        raise ValueError(f"{path.name} must contain sample_id and replicate")
    wave_cols=[]; waves=[]
    for col in df.columns:
        w=parse_wavelength_column(str(col))
        if w is not None: wave_cols.append(col); waves.append(w)
    if len(wave_cols)<100: raise ValueError(f"{path.name}: too few wavelength columns")
    waves_arr=np.asarray(waves,float)
    if not np.all(np.diff(waves_arr)>0): raise ValueError(f"{path.name}: wavelength columns must increase")
    if abs(waves_arr[0]-400.0)>1e-9 or abs(waves_arr[-1]-2499.5)>1e-9:
        raise ValueError(f"{path.name}: expected exported range 400.0-2499.5 nm")
    if abs(float(np.median(np.diff(waves_arr)))-0.5)>1e-9:
        raise ValueError(f"{path.name}: expected exported spacing 0.5 nm")
    df["sample_id"]=df["sample_id"].astype(str)
    replicate=pd.to_numeric(df["replicate"],errors="raise")
    if not replicate.eq(1).all() or df["sample_id"].duplicated().any():
        raise ValueError(f"{path.name}: exactly one observation (replicate=1) per independent sample is required")
    if df["sample_id"].str.strip().isin(['','nan']).any(): raise ValueError('Blank sample ID')
    df=df.sort_values('sample_id').reset_index(drop=True)
    ids=df['sample_id'].to_numpy(str)
    vals=df[wave_cols].to_numpy(dtype=float)
    if not np.isfinite(vals).all(): raise ValueError(f"{path.name}: non-finite spectral value")
    X=vals[:,None,:]
    return ids,waves_arr,X


def load_data(data_dir: Path, grid_step_nm: float, spectral_min_nm: float=400.0, spectral_max_nm: float=2498.0) -> PairedSpectralData:
    whole_path=data_dir/"whole_grain_replicate_spectra.csv.gz"
    ground_path=data_dir/"ground_grain_replicate_spectra.csv.gz"
    chem_path=data_dir/"reference_chemistry.csv"
    for p in (whole_path,ground_path,chem_path):
        if not p.exists(): raise FileNotFoundError(p)
    ids_w,waves_w,Xw=_load_presentation_csv(whole_path)
    ids_g,waves_g,Xg=_load_presentation_csv(ground_path)
    if not np.array_equal(ids_w,ids_g): raise ValueError("Whole/ground sample IDs differ")
    if not np.array_equal(waves_w,waves_g): raise ValueError("Whole/ground wavelength axes differ")
    chem=pd.read_csv(chem_path); expected={"sample_id"}.union(TARGETS)
    if not expected.issubset(chem.columns): raise ValueError(f"reference_chemistry.csv must contain {sorted(expected)}")
    chem["sample_id"]=chem["sample_id"].astype(str)
    if chem["sample_id"].duplicated().any(): raise ValueError("Duplicate sample IDs in reference chemistry")
    if set(chem["sample_id"]) != set(ids_w): raise ValueError("Chemistry/spectral sample IDs differ")
    if not np.isfinite(chem[TARGETS].to_numpy(float)).all(): raise ValueError("Non-finite chemistry")
    chem=chem.set_index("sample_id").loc[ids_w].reset_index()
    if chem[TARGETS].isna().any().any(): raise ValueError("Missing reference chemistry values")
    ratio=grid_step_nm/0.5
    if abs(ratio-round(ratio))>1e-9 or ratio<1: raise ValueError("grid_step_nm must be an integer multiple of 0.5 nm")
    stride=int(round(ratio)); idx=np.arange(0,len(waves_w),stride)
    if spectral_min_nm<400.0-1e-12 or spectral_max_nm>2499.5+1e-12 or spectral_min_nm>=spectral_max_nm:
        raise ValueError("Invalid spectral range")
    idx=idx[(waves_w[idx]>=spectral_min_nm-1e-12)&(waves_w[idx]<=spectral_max_nm+1e-12)]
    if len(idx)<20: raise ValueError("Selected spectral range/grid leaves too few wavelength points")
    return PairedSpectralData(ids_w,waves_w[idx],Xw[:,:,idx],Xg[:,:,idx],chem)


def make_regression_stratified_splits(y: np.ndarray, n_splits: int, n_repeats: int, seed: int,
                                       groups: Optional[np.ndarray]=None) -> List[Tuple[int,int,np.ndarray,np.ndarray]]:
    y=np.asarray(y,float); n=len(y)
    if groups is None: groups=np.asarray([f"sample_{i}" for i in range(n)],dtype=str)
    else:
        groups=np.asarray(groups,dtype=str)
        if len(groups)!=n: raise ValueError("groups must have one value per sample")
    unique_groups=np.array(sorted(np.unique(groups)),dtype=str)
    if n_splits<2 or n_splits>len(unique_groups): raise ValueError("Invalid number of splits for independent CV groups")
    group_y=np.array([float(np.mean(y[groups==g])) for g in unique_groups])
    ordered_groups=unique_groups[np.argsort(group_y,kind="mergesort")]
    out=[]
    for repeat in range(n_repeats):
        rng=np.random.default_rng(stable_seed(seed,"split",repeat,n_splits,n,len(unique_groups)))
        group_fold={}
        for start in range(0,len(ordered_groups),n_splits):
            block=ordered_groups[start:start+n_splits]; labels=np.arange(len(block),dtype=int); rng.shuffle(labels)
            for g,fold in zip(block,labels): group_fold[str(g)]=int(fold)
        fold_id=np.array([group_fold[str(g)] for g in groups],dtype=int)
        for fold in range(n_splits):
            test=np.where(fold_id==fold)[0]; train=np.where(fold_id!=fold)[0]
            if len(test)==0 or len(train)==0: raise RuntimeError("Grouped split produced empty train/test fold")
            out.append((repeat,fold,train,test))
    return out


def load_metadata_groups(data_dir: Path, sample_ids: np.ndarray) -> Tuple[Optional[np.ndarray],Dict[str,Any]]:
    path=data_dir/"sample_metadata_template.csv"
    report={"metadata_file":path.name,"metadata_file_present":path.exists(),"source_group_status":"missing_file",
            "cv_grouping_used":"sample_id","submission_warning":None}
    if not path.exists():
        report["submission_warning"]="Sample provenance metadata is absent; independence beyond coded sample IDs is unverified."
        return None,report
    df=pd.read_csv(path,dtype=str).fillna("")
    if "sample_id" not in df.columns: raise ValueError("sample_metadata_template.csv must contain sample_id")
    df["sample_id"]=df["sample_id"].astype(str)
    if df["sample_id"].duplicated().any(): raise ValueError("Duplicate sample IDs in sample metadata")
    if set(df["sample_id"])!=set(sample_ids.astype(str)): raise ValueError("Sample metadata IDs do not match spectral IDs")
    df=df.set_index("sample_id").loc[sample_ids.astype(str)].reset_index()
    report["field_completeness_fraction"]={c:float(np.mean(df[c].astype(str).str.strip()!="")) for c in df.columns if c!="sample_id"}
    if "source_group" not in df.columns:
        report["source_group_status"]="column_absent"; report["submission_warning"]="source_group is absent; coded samples are treated as independent."
        return None,report
    values=df["source_group"].astype(str).str.strip(); n_filled=int(np.sum(values!=""))
    if n_filled==0:
        report["source_group_status"]="unverified_blank"
        report["submission_warning"]="source_group is blank. CV treats 98 user-declared independent samples as independent; complete provenance before submission if any samples share field/batch/material origin."
        return None,report
    if n_filled!=len(values): raise ValueError("source_group is partially populated; complete all rows or leave all blank")
    groups=values.to_numpy(str); report["source_group_status"]="complete"; report["n_unique_source_groups"]=int(len(np.unique(groups)))
    report["cv_grouping_used"]="source_group" if len(np.unique(groups))<len(groups) else "sample_id_equivalent_source_group"
    return groups,report
