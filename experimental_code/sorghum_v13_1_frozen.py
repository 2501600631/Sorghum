from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
sys.dont_write_bytecode = True
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.signal import savgol_filter
from scipy.linalg import qr
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR


from sorghum_v13_data_utils import (
    TARGETS,
    TARGET_LABELS,
    PairedSpectralData,
    ensure_dir,
    load_data,
    load_metadata_groups,
    make_regression_stratified_splits,
    regression_metrics,
    save_json,
    sha256_file,
    stable_seed,
)

SCRIPT_VERSION = "13.2.0-98-final-methodological-hardening"
DEFAULT_SEED = 20260830
FRAMEWORK_VERSION = '13.2.0-98'
PRIMARY_TARGET = "crude_protein"
SECONDARY_TARGET = "moisture"






@dataclass(frozen=True)
class Representation:
    name: str
    min_nm: float
    max_nm: float
    preprocess: str
    sg_window: int = 17
    sg_poly: int = 2


def _snv(X: np.ndarray) -> np.ndarray:
    mu = np.mean(X, axis=1, keepdims=True)
    sd = np.std(X, axis=1, ddof=1, keepdims=True)
    sd = np.where(sd < 1e-12, 1.0, sd)
    return (X - mu) / sd


class SpectralPreprocessor:





    def __init__(self, method: str, wavelengths: np.ndarray, sg_window: int = 17, sg_poly: int = 2):
        self.method = str(method).upper()
        self.wavelengths = np.asarray(wavelengths, float)
        self.sg_window = int(sg_window)
        self.sg_poly = int(sg_poly)
        self.msc_reference_: Optional[np.ndarray] = None

    def _valid_window(self, p: int) -> int:
        w = min(self.sg_window, p if p % 2 == 1 else p - 1)
        w = max(w, self.sg_poly + 2 + ((self.sg_poly + 2) % 2 == 0))
        if w % 2 == 0:
            w += 1
        if w > p:
            w = p if p % 2 == 1 else p - 1
        return max(w, 5)

    def _base_transform(self, X: np.ndarray, token: str) -> np.ndarray:
        X = np.asarray(X, float)
        token = token.upper()
        if token == "RAW":
            return X.copy()
        if token == "SNV":
            return _snv(X)
        if token == "SG":
            w = self._valid_window(X.shape[1])
            return savgol_filter(X, window_length=w, polyorder=min(self.sg_poly, w - 2), deriv=0, axis=1, mode="interp")
        if token == "SG1":
            w = self._valid_window(X.shape[1])
            delta = float(np.median(np.diff(self.wavelengths)))
            return savgol_filter(X, window_length=w, polyorder=min(max(self.sg_poly, 2), w - 2), deriv=1, delta=delta, axis=1, mode="interp")
        if token == "SG2":
            w = self._valid_window(X.shape[1])
            delta = float(np.median(np.diff(self.wavelengths)))
            return savgol_filter(X, window_length=w, polyorder=min(max(self.sg_poly, 3), w - 2), deriv=2, delta=delta, axis=1, mode="interp")
        if token == "DE1":
            return np.gradient(X, self.wavelengths, axis=1)
        if token == "DE2":
            return np.gradient(np.gradient(X, self.wavelengths, axis=1), self.wavelengths, axis=1)
        raise ValueError(f"Unknown preprocessing token: {token}")

    @staticmethod
    def _msc_apply(X: np.ndarray, ref: np.ndarray) -> np.ndarray:
        ref = np.asarray(ref, float)
        rc = ref - ref.mean()
        den = float(np.dot(rc, rc))
        if den < 1e-20:
            return X.copy()
        out = np.empty_like(X, dtype=float)
        rmean = float(ref.mean())
        for i, row in enumerate(np.asarray(X, float)):
            slope = float(np.dot(row - row.mean(), rc) / den)
            if abs(slope) < 1e-12:
                slope = 1.0
            intercept = float(row.mean() - slope * rmean)
            out[i] = (row - intercept) / slope
        return out

    def _before_msc(self, X: np.ndarray) -> np.ndarray:
        m = self.method
        if m == "MSC":
            return np.asarray(X, float)
        if m.endswith("-MSC"):
            first = m[:-4]
            return self._base_transform(X, first)
        raise ValueError("MSC requested for non-MSC method")

    def fit(self, X: np.ndarray) -> "SpectralPreprocessor":
        if "MSC" in self.method:
            base = self._before_msc(X)
            self.msc_reference_ = np.mean(base, axis=0)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        m = self.method
        if m == "RAW":
            return np.asarray(X, float).copy()
        if m in {"SNV", "SG", "SG1", "SG2", "DE1", "DE2"}:
            return self._base_transform(X, m)
        if m == "MSC" or m.endswith("-MSC"):
            if self.msc_reference_ is None:
                raise RuntimeError("MSC preprocessor must be fitted before transform")
            base = self._before_msc(X)
            return self._msc_apply(base, self.msc_reference_)
        if m.endswith("-SNV"):
            first = m[:-4]
            return _snv(self._base_transform(X, first))
        raise ValueError(f"Unsupported preprocessing method: {m}")


def representation_space(profile: str) -> List[Representation]:
    profile = profile.lower()
    if profile == "smoke":

        ranges = [("FULL", 400.0, 2498.0)]
        methods = ["SNV"]
        windows = [(17, 2)]
    elif profile == "extended":
        ranges = [
            ("FULL", 400.0, 2498.0), ("NIR900", 900.0, 2498.0),
            ("NIR1100", 1100.0, 2498.0), ("NIR1100_1800", 1100.0, 1800.0),
            ("NIR1800", 1800.0, 2498.0), ("PROTEIN_A", 900.0, 1300.0),
        ]
        methods = ["RAW", "SNV", "MSC", "SG", "SG1-SNV", "SG1-MSC", "DE1-SNV", "DE1-MSC", "DE2-SNV"]
        windows = [(11, 2), (17, 2), (21, 2), (21, 3)]
    else:
        ranges = [
            ("FULL", 400.0, 2498.0), ("NIR900", 900.0, 2498.0),
            ("NIR1100", 1100.0, 2498.0), ("NIR1100_1800", 1100.0, 1800.0),
            ("NIR1800", 1800.0, 2498.0),
        ]

        methods = ["RAW", "SNV", "MSC", "SG1-SNV", "DE1-SNV", "DE1-MSC"]
        windows = [(17, 2)]
    reps: List[Representation] = []
    for rn, lo, hi in ranges:
        for method in methods:
            wset = windows if ("SG" in method or profile == "extended") else [(17, 2)]
            for win, poly in wset:
                suffix = f"_w{win}p{poly}" if "SG" in method else ""
                reps.append(Representation(f"{rn}__{method}{suffix}", lo, hi, method, win, poly))

    uniq: Dict[str, Representation] = {}
    for r in reps:
        uniq[r.name] = r
    return list(uniq.values())


def rep_mask(wavelengths: np.ndarray, rep: Representation) -> np.ndarray:
    return (wavelengths >= rep.min_nm - 1e-12) & (wavelengths <= rep.max_nm + 1e-12)


def transform_rep(rep: Representation, wavelengths: np.ndarray, Xw_fit: np.ndarray, *arrays: np.ndarray) -> Tuple[SpectralPreprocessor, List[np.ndarray]]:
    mask = rep_mask(wavelengths, rep)
    if int(mask.sum()) < 20:
        raise ValueError(f"Representation {rep.name} leaves too few wavelengths")
    pp = SpectralPreprocessor(rep.preprocess, wavelengths[mask], rep.sg_window, rep.sg_poly).fit(Xw_fit[:, mask])
    out = [pp.transform(a[:, mask]) for a in arrays]
    return pp, out






def _standardize_fit(X: np.ndarray) -> Tuple[StandardScaler, np.ndarray]:
    sc = StandardScaler().fit(X)
    return sc, sc.transform(X)


def _fit_pls(X: np.ndarray, y: np.ndarray, n_components: int) -> Dict[str, Any]:
    ncomp = max(1, min(int(n_components), X.shape[0] - 1, X.shape[1]))
    m = PLSRegression(n_components=ncomp, scale=True).fit(X, y)
    return {"model": m}


def _fit_ridge(X: np.ndarray, y: np.ndarray, alpha: float) -> Dict[str, Any]:
    sc, Z = _standardize_fit(X)
    m = Ridge(alpha=float(alpha)).fit(Z, y)
    return {"scaler": sc, "model": m}


class ElasticNetConvergenceError(RuntimeError):

    def __init__(self, diagnostic: Dict[str, Any]):
        self.diagnostic = diagnostic
        super().__init__("ElasticNet candidate failed the convergence-safe fitting policy")


def _fit_elastic(X: np.ndarray, y: np.ndarray, alpha_ratio: float, l1_ratio: float, seed: int) -> Dict[str, Any]:







    sc, Z = _standardize_fit(X)
    yy = np.asarray(y, float).reshape(-1)
    y_mean = float(np.mean(yy))
    y_scale = float(np.std(yy, ddof=0))
    if not np.isfinite(y_scale) or y_scale < 1e-12:
        y_scale = 1.0
    yz = (yy - y_mean) / y_scale
    l1 = float(l1_ratio)
    if not (0.0 < l1 <= 1.0):
        raise ValueError("ElasticNet l1_ratio must lie in (0, 1]")
    alpha_max = float(np.max(np.abs(Z.T @ yz)) / (max(len(yz), 1) * l1))
    alpha_max = max(alpha_max, np.finfo(float).eps)
    ratio = float(alpha_ratio)
    if not (0.0 < ratio <= 1.0):
        raise ValueError("ElasticNet alpha_ratio must lie in (0, 1]")
    alpha = ratio * alpha_max

    attempts = [
        {"max_iter": 100000, "tol": 1e-4, "selection": "random"},
        {"max_iter": 500000, "tol": 5e-5, "selection": "random"},
    ]
    attempt_records: List[Dict[str, Any]] = []
    for attempt_id, cfg in enumerate(attempts, start=1):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ConvergenceWarning)
            model = ElasticNet(
                alpha=alpha, l1_ratio=l1, max_iter=int(cfg["max_iter"]),
                tol=float(cfg["tol"]), selection=str(cfg["selection"]),
                random_state=stable_seed(seed, "elasticnet_solver", attempt_id),
                precompute=True,
            ).fit(Z, yz)
        had_warning = any(issubclass(w.category, ConvergenceWarning) for w in caught)
        n_iter = int(np.max(np.atleast_1d(model.n_iter_)))
        dual_gap = float(np.max(np.atleast_1d(model.dual_gap_)))
        converged = bool((not had_warning) and np.isfinite(dual_gap) and n_iter < int(cfg["max_iter"]))
        attempt_records.append({
            "attempt": attempt_id, "max_iter": int(cfg["max_iter"]), "tol": float(cfg["tol"]),
            "selection": str(cfg["selection"]), "n_iter": n_iter, "dual_gap": dual_gap,
            "convergence_warning": bool(had_warning), "converged": converged,
        })
        if converged:
            return {
                "scaler": sc, "model": model, "y_mean": y_mean, "y_scale": y_scale,
                "elasticnet_diagnostic": {
                    "converged": True, "alpha_ratio": ratio, "alpha_max": alpha_max,
                    "alpha_absolute": alpha, "l1_ratio": l1, "attempts_used": attempt_id,
                    "final_n_iter": n_iter, "final_dual_gap": dual_gap,
                    "attempt_records": attempt_records,
                },
            }
    diagnostic = {
        "converged": False, "alpha_ratio": ratio, "alpha_max": alpha_max,
        "alpha_absolute": alpha, "l1_ratio": l1, "attempts_used": len(attempts),
        "final_n_iter": attempt_records[-1]["n_iter"],
        "final_dual_gap": attempt_records[-1]["dual_gap"],
        "attempt_records": attempt_records,
    }
    raise ElasticNetConvergenceError(diagnostic)


def _fit_svr(X: np.ndarray, y: np.ndarray, C: float, gamma: float, epsilon: float) -> Dict[str, Any]:


    sc, Z = _standardize_fit(X)
    yy = np.asarray(y, float).reshape(-1)
    y_mean = float(np.mean(yy)); y_scale = float(np.std(yy, ddof=0))
    if not np.isfinite(y_scale) or y_scale < 1e-12:
        y_scale = 1.0
    yz = (yy - y_mean) / y_scale
    m = SVR(kernel="rbf", C=float(C), gamma=float(gamma), epsilon=float(epsilon)).fit(Z, yz)
    return {"scaler": sc, "model": m, "y_mean": y_mean, "y_scale": y_scale}


def _fit_krr(X: np.ndarray, y: np.ndarray, alpha: float, gamma: float) -> Dict[str, Any]:
    sc, Z = _standardize_fit(X)
    m = KernelRidge(alpha=float(alpha), kernel="rbf", gamma=float(gamma)).fit(Z, y)
    return {"scaler": sc, "model": m}


def _fit_extratrees(X: np.ndarray, y: np.ndarray, max_features: float, min_samples_leaf: int, seed: int) -> Dict[str, Any]:
    m = ExtraTreesRegressor(
        n_estimators=400, max_features=float(max_features), min_samples_leaf=int(min_samples_leaf),
        random_state=int(seed), n_jobs=1,
    ).fit(X, y)
    return {"model": m}


def _fit_pls_krr(X: np.ndarray, y: np.ndarray, n_components: int, alpha: float, gamma: float) -> Dict[str, Any]:
    ncomp = max(1, min(int(n_components), X.shape[0] - 1, X.shape[1]))
    pls = PLSRegression(n_components=ncomp, scale=True).fit(X, y)
    T = pls.x_scores_
    sc = StandardScaler().fit(T)
    krr = KernelRidge(alpha=float(alpha), kernel="rbf", gamma=float(gamma)).fit(sc.transform(T), y)
    return {"pls": pls, "score_scaler": sc, "model": krr}


def _predict_basic(kind: str, obj: Dict[str, Any], X: np.ndarray) -> np.ndarray:
    if kind in {"pls", "extratrees"}:
        return np.asarray(obj["model"].predict(X)).reshape(-1)
    if kind == "elastic":
        z = np.asarray(obj["model"].predict(obj["scaler"].transform(X))).reshape(-1)
        return obj["y_mean"] + obj["y_scale"] * z
    if kind == "svr":
        z = np.asarray(obj["model"].predict(obj["scaler"].transform(X))).reshape(-1)
        return obj["y_mean"] + obj["y_scale"] * z
    if kind in {"ridge", "krr"}:
        return np.asarray(obj["model"].predict(obj["scaler"].transform(X))).reshape(-1)
    if kind == "pls_krr":
        T = obj["pls"].transform(X)
        return np.asarray(obj["model"].predict(obj["score_scaler"].transform(T))).reshape(-1)
    raise ValueError(kind)


def _fit_basic(kind: str, X: np.ndarray, y: np.ndarray, params: Dict[str, Any], seed: int) -> Dict[str, Any]:
    if kind == "pls": return _fit_pls(X, y, params["n_components"])
    if kind == "ridge": return _fit_ridge(X, y, params["alpha"])
    if kind == "elastic": return _fit_elastic(X, y, params["alpha_ratio"], params["l1_ratio"], seed)
    if kind == "svr": return _fit_svr(X, y, params["C"], params["gamma"], params["epsilon"])
    if kind == "krr": return _fit_krr(X, y, params["alpha"], params["gamma"])
    if kind == "extratrees": return _fit_extratrees(X, y, params["max_features"], params["min_samples_leaf"], seed)
    if kind == "pls_krr": return _fit_pls_krr(X, y, params["n_components"], params["alpha"], params["gamma"])
    raise ValueError(kind)


class IPLSModel:

    def __init__(self, n_intervals: int, n_components: int, seed: int):
        self.n_intervals = int(n_intervals)
        self.n_components = int(n_components)
        self.seed = int(seed)
        self.interval_: Optional[np.ndarray] = None
        self.model_: Optional[PLSRegression] = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> "IPLSModel":
        p = X.shape[1]
        intervals = [idx for idx in np.array_split(np.arange(p), self.n_intervals) if len(idx) >= 3]
        splits = make_regression_stratified_splits(y, min(3, max(2, len(y)//8)), 1, self.seed)
        best_rmse = float("inf"); best_idx = None
        for idx in intervals:
            pred = np.full(len(y), np.nan)
            ok = True
            for _, _, tr, va in splits:
                nc = min(self.n_components, len(tr)-1, len(idx))
                if nc < 1:
                    ok = False; break
                try:
                    m = PLSRegression(n_components=nc, scale=True).fit(X[tr][:, idx], y[tr])
                    pred[va] = m.predict(X[va][:, idx]).reshape(-1)
                except (np.linalg.LinAlgError, FloatingPointError):
                    ok = False; break
            if ok and np.all(np.isfinite(pred)):
                rmse = float(np.sqrt(np.mean((pred-y)**2)))
                if rmse < best_rmse:
                    best_rmse = rmse; best_idx = idx
        if best_idx is None:
            raise RuntimeError("iPLS could not fit any interval")
        nc = min(self.n_components, len(y)-1, len(best_idx))
        self.interval_ = np.asarray(best_idx, int)
        self.model_ = PLSRegression(n_components=nc, scale=True).fit(X[:, self.interval_], y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if self.interval_ is None or self.model_ is None:
            raise RuntimeError("iPLS model not fitted")
        return self.model_.predict(X[:, self.interval_]).reshape(-1)







def fit_epo_base(Xw: np.ndarray, Xg: np.ndarray, y: np.ndarray, n_epo: int, base_kind: str, base_params: Dict[str, Any], seed: int) -> Dict[str, Any]:
    D = np.asarray(Xg - Xw, float)
    if n_epo <= 0:
        V = np.empty((Xw.shape[1], 0))
    else:
        _, _, vt = np.linalg.svd(D, full_matrices=False)
        k = min(int(n_epo), vt.shape[0])
        V = vt[:k].T
    def proj(X):
        return X if V.shape[1] == 0 else X - (X @ V) @ V.T
    base = _fit_basic(base_kind, proj(Xw), y, base_params, seed)
    return {"V": V, "base_kind": base_kind, "base": base}


def predict_epo(obj: Dict[str, Any], Xw: np.ndarray) -> np.ndarray:
    V = obj["V"]
    Z = Xw if V.shape[1] == 0 else Xw - (Xw @ V) @ V.T
    return _predict_basic(obj["base_kind"], obj["base"], Z)



def fit_pooled_base(Xw: np.ndarray, Xg: np.ndarray, y: np.ndarray, base_kind: str, base_params: Dict[str, Any], seed: int) -> Dict[str, Any]:
    X = np.vstack([Xw, Xg]); yy = np.concatenate([y, y])
    return {"base_kind": base_kind, "base": _fit_basic(base_kind, X, yy, base_params, seed)}


def predict_pooled(obj: Dict[str, Any], Xw: np.ndarray) -> np.ndarray:
    return _predict_basic(obj["base_kind"], obj["base"], Xw)



def _domain_penalty_basis(xs: np.ndarray, xt: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    xs = np.asarray(xs, float); xt = np.asarray(xt, float)
    xs = xs - xs.mean(axis=0, keepdims=True)
    xt = xt - xt.mean(axis=0, keepdims=True)
    A = np.vstack([xs / math.sqrt(max(len(xs),1)), xt / math.sqrt(max(len(xt),1))])
    signs = np.concatenate([np.ones(len(xs)), -np.ones(len(xt))])

    Q, R = qr(A.T, mode="economic")
    M = (R * signs[np.newaxis, :]) @ R.T
    evals, U = np.linalg.eigh(0.5*(M+M.T))
    keep = np.abs(evals) > 1e-12
    if not np.any(keep):
        return np.empty((xs.shape[1],0)), np.empty(0)
    V = Q @ U[:, keep]
    s = np.abs(evals[keep])
    return V, s


def fit_dipls_whole_target(Xw: np.ndarray, Xg: np.ndarray, y: np.ndarray, n_components: int, lam: float) -> Dict[str, Any]:






    Xmean = Xw.mean(axis=0); Gmean = Xg.mean(axis=0); ymean = float(np.mean(y))
    x = Xw - Xmean; xs = x.copy(); xt = Xg - Gmean; yy = (y-ymean).reshape(-1,1)
    p = x.shape[1]; A = min(int(n_components), len(y)-1, p)
    W=[]; P=[]; C=[]
    Ieps = 1e-12
    for _ in range(A):
        denom_y = float((yy.T @ yy).item())
        if denom_y < 1e-18: break
        w_pls = (yy.T @ x / denom_y).reshape(-1)
        norm = np.linalg.norm(w_pls)
        if norm < 1e-18: break
        if float(lam) > 0:
            V,s = _domain_penalty_basis(xs, xt)
            w = w_pls.copy()
            if len(s):
                c = float(lam)/denom_y
                coeff = (c*s)/(1.0+c*s)
                w = w - V @ (coeff * (V.T @ w))
        else:
            w = w_pls.copy()
        nw=np.linalg.norm(w)
        if nw < 1e-18: break
        w=w/nw
        t=x@w.reshape(-1,1); ts=xs@w.reshape(-1,1); tt=xt@w.reshape(-1,1)
        tden=float((t.T@t).item()); tsden=float((ts.T@ts).item()); ttden=float((tt.T@tt).item())
        if tden<Ieps or tsden<Ieps or ttden<Ieps: break
        ccoef=float((yy.T@t).item()/tden)
        pvec=(t.T@x/tden).reshape(-1)
        ps=(ts.T@xs/tsden).reshape(-1); pt=(tt.T@xt/ttden).reshape(-1)
        x=x-t@pvec.reshape(1,-1); xs=xs-ts@ps.reshape(1,-1); xt=xt-tt@pt.reshape(1,-1); yy=yy-t*ccoef
        W.append(w); P.append(pvec); C.append(ccoef)
    if not W:
        raise RuntimeError("di-PLS produced no latent variables")
    Wm=np.column_stack(W); Pm=np.column_stack(P); Cm=np.asarray(C).reshape(-1,1)
    b=Wm@np.linalg.pinv(Pm.T@Wm)@Cm
    return {"coef": b.reshape(-1), "x_mean": Xmean, "y_mean": ymean, "n_components_fitted": len(W)}


def predict_dipls(obj: Dict[str, Any], Xw: np.ndarray) -> np.ndarray:
    return obj["y_mean"] + (Xw-obj["x_mean"])@obj["coef"]






@dataclass(frozen=True)
class Candidate:
    family: str
    rep_name: str
    params_json: str
    complexity: int

    @property
    def params(self) -> Dict[str, Any]:
        return json.loads(self.params_json)

    @property
    def key(self) -> str:
        h = hashlib.sha1((self.family+"|"+self.rep_name+"|"+self.params_json).encode()).hexdigest()[:10]
        return f"{self.family}__{self.rep_name}__{h}"


def C(family: str, rep: Representation, params: Dict[str, Any], complexity: int) -> Candidate:
    return Candidate(family, rep.name, json.dumps(params, sort_keys=True, separators=(",",":")), complexity)


def _candidate_grids(rep: Representation, profile: str, transfer: bool) -> List[Candidate]:
    smoke = profile == "smoke"
    out: List[Candidate]=[]
    if not transfer:
        pls=[2] if smoke else [2,3,4,5,6,8,10,12,15]
        ridge=[1] if smoke else [1e-4,1e-3,0.01,0.1,1,10,100,1000]
        en_a=[0.1] if smoke else [0.3,0.1,0.03,0.01,0.003]
        en_l=[0.5] if smoke else [0.1,0.5,0.9]
        for n in pls: out.append(C("whole_pls",rep,{"n_components":n},1))
        for a in ridge: out.append(C("whole_ridge",rep,{"alpha":a},1))
        for a in en_a:
            for l in en_l: out.append(C("whole_elastic",rep,{"alpha_ratio":a,"l1_ratio":l},2))
        if not smoke:
            for Cc in [1,10,100,1000]:
                for g in [1e-3,1e-2,0.1]:
                    for e in [0.05,0.1,0.2]: out.append(C("whole_svr",rep,{"C":Cc,"gamma":g,"epsilon":e},4))
            for a in [0.01,0.1,1,10]:
                for g in [1e-4,1e-3,1e-2,0.1]: out.append(C("whole_krr",rep,{"alpha":a,"gamma":g},4))
            for mf in [0.3,0.6,1.0]:
                for leaf in [1,2]: out.append(C("whole_extratrees",rep,{"max_features":mf,"min_samples_leaf":leaf},5))
            for n in [4,6,8,10]:
                for a in [0.01,0.1,1]:
                    for g in [0.01,0.1]: out.append(C("whole_pls_krr",rep,{"n_components":n,"alpha":a,"gamma":g},5))
            for ni in [10,20,30]:
                for n in [4,6,8,10]: out.append(C("whole_ipls",rep,{"n_intervals":ni,"n_components":n},3))
    else:
        pls=[2] if smoke else [2,4,6,8,10]
        epo=[1] if smoke else [1,2,3,5,8]
        for ne in epo:
            for n in pls: out.append(C("epo_pls",rep,{"n_epo":ne,"n_components":n},3))
        en_a=[0.1] if smoke else [0.1,0.03,0.01,0.003]
        en_l=[0.5] if smoke else [0.1,0.5,0.9]
        for ne in epo:
            for a in en_a:
                for l in en_l: out.append(C("epo_elastic",rep,{"n_epo":ne,"alpha_ratio":a,"l1_ratio":l},4))
        for a in en_a:
            for l in en_l: out.append(C("pooled_elastic",rep,{"alpha_ratio":a,"l1_ratio":l},3))
        dip_n=[2] if smoke else [2,4,6,8]
        dip_l=[1.0] if smoke else [0.0,0.1,1,10,100,1000,10000,100000]
        for n in dip_n:
            for lam in dip_l: out.append(C("dipls",rep,{"n_components":n,"lambda":lam},4))
    return out


def _family_base_kind(family: str) -> str:
    return {
        "whole_pls":"pls","whole_ridge":"ridge","whole_elastic":"elastic","whole_svr":"svr",
        "whole_krr":"krr","whole_extratrees":"extratrees","whole_pls_krr":"pls_krr"
    }[family]


def fit_candidate(cand: Candidate, rep: Representation, wavelengths: np.ndarray,
                  Xw_train_raw: np.ndarray, Xg_train_raw: np.ndarray, y: np.ndarray, seed: int) -> Dict[str, Any]:
    mask=rep_mask(wavelengths,rep)
    pp=SpectralPreprocessor(rep.preprocess,wavelengths[mask],rep.sg_window,rep.sg_poly).fit(Xw_train_raw[:,mask])
    Xw=pp.transform(Xw_train_raw[:,mask]); Xg=pp.transform(Xg_train_raw[:,mask])
    p=cand.params; f=cand.family
    if f.startswith("whole_") and f!="whole_ipls":
        kind=_family_base_kind(f); obj=_fit_basic(kind,Xw,y,p,seed)
        return {"rep":rep,"mask":mask,"pp":pp,"family":f,"kind":kind,"obj":obj}
    if f=="whole_ipls":
        m=IPLSModel(p["n_intervals"],p["n_components"],stable_seed(seed,"ipls")).fit(Xw,y)
        return {"rep":rep,"mask":mask,"pp":pp,"family":f,"obj":m}
    if f=="epo_pls":
        obj=fit_epo_base(Xw,Xg,y,p["n_epo"],"pls",{"n_components":p["n_components"]},seed)
    elif f=="epo_elastic":
        obj=fit_epo_base(Xw,Xg,y,p["n_epo"],"elastic",{"alpha_ratio":p["alpha_ratio"],"l1_ratio":p["l1_ratio"]},seed)
    elif f=="pooled_elastic":
        obj=fit_pooled_base(Xw,Xg,y,"elastic",{"alpha_ratio":p["alpha_ratio"],"l1_ratio":p["l1_ratio"]},seed)
    elif f=="dipls":
        obj=fit_dipls_whole_target(Xw,Xg,y,p["n_components"],p["lambda"])
    else:
        raise ValueError(f)
    return {"rep":rep,"mask":mask,"pp":pp,"family":f,"obj":obj}


def predict_candidate(fit: Dict[str,Any], Xw_raw: np.ndarray) -> np.ndarray:
    mask=fit["mask"]; Xw=fit["pp"].transform(Xw_raw[:,mask]); f=fit["family"]
    if f.startswith("whole_") and f!="whole_ipls": return _predict_basic(fit["kind"],fit["obj"],Xw)
    if f=="whole_ipls": return fit["obj"].predict(Xw)
    if f in {"epo_pls","epo_elastic"}: return predict_epo(fit["obj"],Xw)
    if f=="pooled_elastic": return predict_pooled(fit["obj"],Xw)
    if f=="dipls": return predict_dipls(fit["obj"],Xw)
    raise ValueError(f)


def _extract_elasticnet_diagnostic(fit: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    f = fit.get("family", "")
    if f == "whole_elastic":
        return fit["obj"].get("elasticnet_diagnostic")
    if f in {"epo_elastic", "pooled_elastic"}:
        return fit["obj"]["base"].get("elasticnet_diagnostic")
    return None


def cv_candidate(cand: Candidate, rep_lookup: Dict[str,Representation], wavelengths: np.ndarray,
                 Xw: np.ndarray, Xg: np.ndarray, y: np.ndarray,
                 inner: List[Tuple[np.ndarray,np.ndarray]], seed: int) -> Tuple[float,float,List[float],List[Dict[str,Any]],str,str,str]:







    rmses=[]; elastic_audit: List[Dict[str,Any]]=[]
    for j,(tr,va) in enumerate(inner):
        try:
            fit=fit_candidate(cand,rep_lookup[cand.rep_name],wavelengths,Xw[tr],Xg[tr],y[tr],stable_seed(seed,cand.key,j))
            diag=_extract_elasticnet_diagnostic(fit)
            if diag is not None:
                d=dict(diag); d["inner_fold"]=int(j); d["candidate_key"]=cand.key; d["family"]=cand.family; elastic_audit.append(d)
            pred=predict_candidate(fit,Xw[va])
            if not np.all(np.isfinite(pred)):
                return float("inf"),float("inf"),[],elastic_audit,"invalid","NonFinitePrediction",f"non-finite prediction in inner fold {j}"
            rmses.append(float(np.sqrt(np.mean((pred-y[va])**2))))
        except ElasticNetConvergenceError as exc:
            d=dict(exc.diagnostic); d["inner_fold"]=int(j); d["candidate_key"]=cand.key; d["family"]=cand.family; elastic_audit.append(d)
            return float("inf"),float("inf"),[],elastic_audit,"invalid","ElasticNetConvergenceError",str(exc)
        except (np.linalg.LinAlgError, FloatingPointError) as exc:
            return float("inf"),float("inf"),[],elastic_audit,"invalid",type(exc).__name__,str(exc)
    if len(rmses)!=len(inner):
        return float("inf"),float("inf"),[],elastic_audit,"invalid","IncompleteInnerCV","not all inner folds produced predictions"
    return float(np.mean(rmses)),float(np.std(rmses,ddof=1)/math.sqrt(len(rmses))) if len(rmses)>1 else 0.0,rmses,elastic_audit,"valid","",""

def select_min_inner_rmse(records: pd.DataFrame) -> pd.Series:







    good=records[np.isfinite(records["mean_rmse"])].copy()
    if good.empty: raise RuntimeError("No valid candidate in inner CV")
    ordered=good.sort_values(["mean_rmse","candidate_key"], kind="mergesort")
    return ordered.iloc[0]


def _screen_proxy_candidates(rep: Representation, smoke: bool, branch: str) -> List[Candidate]:














    branch=str(branch).lower()
    if branch not in {"whole","transfer","ground"}:
        raise ValueError(branch)
    if smoke:
        if branch == "whole":
            return [
                C("whole_pls", rep, {"n_components":2}, 1),
                C("whole_elastic", rep, {"alpha_ratio":0.1,"l1_ratio":0.5}, 2),
                C("whole_svr", rep, {"C":10.0,"gamma":0.01,"epsilon":0.1}, 4),
                C("whole_krr", rep, {"alpha":0.1,"gamma":0.01}, 4),
            ]
        if branch == "ground":
            return [
                C("ground_pls", rep, {"n_components":2}, 1),
                C("ground_elastic", rep, {"alpha_ratio":0.1,"l1_ratio":0.5}, 2),
                C("ground_svr", rep, {"C":10.0,"gamma":0.01,"epsilon":0.1}, 4),
                C("ground_krr", rep, {"alpha":0.1,"gamma":0.01}, 4),
            ]
        return [
            C("epo_pls", rep, {"n_epo":1,"n_components":2}, 3),
            C("epo_elastic", rep, {"n_epo":1,"alpha_ratio":0.1,"l1_ratio":0.5}, 4),
            C("pooled_elastic", rep, {"alpha_ratio":0.1,"l1_ratio":0.5}, 3),
            C("dipls", rep, {"n_components":2,"lambda":1.0}, 4),
        ]
    if branch == "whole":
        return [
            C("whole_pls", rep, {"n_components":6}, 1),
            C("whole_elastic", rep, {"alpha_ratio":0.03,"l1_ratio":0.5}, 2),
            C("whole_svr", rep, {"C":100.0,"gamma":0.01,"epsilon":0.1}, 4),
            C("whole_krr", rep, {"alpha":0.1,"gamma":0.01}, 4),
        ]
    if branch == "ground":
        return [
            C("ground_pls", rep, {"n_components":6}, 1),
            C("ground_elastic", rep, {"alpha_ratio":0.03,"l1_ratio":0.5}, 2),
            C("ground_svr", rep, {"C":100.0,"gamma":0.01,"epsilon":0.1}, 4),
            C("ground_krr", rep, {"alpha":0.1,"gamma":0.01}, 4),
        ]
    return [
        C("epo_pls", rep, {"n_epo":2,"n_components":6}, 3),
        C("epo_elastic", rep, {"n_epo":2,"alpha_ratio":0.03,"l1_ratio":0.5}, 4),
        C("pooled_elastic", rep, {"alpha_ratio":0.03,"l1_ratio":0.5}, 3),
        C("dipls", rep, {"n_components":6,"lambda":10.0}, 4),
    ]


def _transfer_screen_candidates(rep: Representation, smoke: bool) -> List[Candidate]:

    return _screen_proxy_candidates(rep, smoke, "transfer")

def _ground_family_name(whole_family: str) -> str:
    if not whole_family.startswith("whole_"):
        raise ValueError(whole_family)
    return "ground_" + whole_family[len("whole_"):]


def _ground_candidate_grids(rep: Representation, profile: str) -> List[Candidate]:

    out=[]
    for c in _candidate_grids(rep, profile, transfer=False):
        out.append(C(_ground_family_name(c.family), rep, c.params, c.complexity))
    return out


def fit_ground_candidate(cand: Candidate, rep: Representation, wavelengths: np.ndarray,
                         Xg_train_raw: np.ndarray, y: np.ndarray, seed: int) -> Dict[str, Any]:

    if not cand.family.startswith("ground_"):
        raise ValueError(cand.family)
    mask=rep_mask(wavelengths,rep)
    pp=SpectralPreprocessor(rep.preprocess,wavelengths[mask],rep.sg_window,rep.sg_poly).fit(Xg_train_raw[:,mask])
    Xg=pp.transform(Xg_train_raw[:,mask]); p=cand.params; suffix=cand.family[len("ground_"):]
    if suffix == "ipls":
        m=IPLSModel(p["n_intervals"],p["n_components"],stable_seed(seed,"ground_ipls")).fit(Xg,y)
        return {"rep":rep,"mask":mask,"pp":pp,"family":cand.family,"obj":m}
    kind={
        "pls":"pls","ridge":"ridge","elastic":"elastic","svr":"svr",
        "krr":"krr","extratrees":"extratrees","pls_krr":"pls_krr",
    }[suffix]
    obj=_fit_basic(kind,Xg,y,p,seed)
    return {"rep":rep,"mask":mask,"pp":pp,"family":cand.family,"kind":kind,"obj":obj}


def predict_ground_candidate(fit: Dict[str,Any], Xg_raw: np.ndarray) -> np.ndarray:
    Xg=fit["pp"].transform(Xg_raw[:,fit["mask"]]); suffix=fit["family"][len("ground_"):]
    if suffix == "ipls": return fit["obj"].predict(Xg)
    return _predict_basic(fit["kind"],fit["obj"],Xg)


def _extract_ground_elasticnet_diagnostic(fit: Dict[str,Any]) -> Optional[Dict[str,Any]]:
    if fit.get("family") == "ground_elastic":
        return fit["obj"].get("elasticnet_diagnostic")
    return None


def cv_ground_candidate(cand: Candidate, rep_lookup: Dict[str,Representation], wavelengths: np.ndarray,
                        Xg: np.ndarray, y: np.ndarray, inner: List[Tuple[np.ndarray,np.ndarray]],
                        seed: int) -> Tuple[float,float,List[float],List[Dict[str,Any]],str,str,str]:
    rmses=[]; elastic_audit=[]
    for j,(tr,va) in enumerate(inner):
        try:
            fit=fit_ground_candidate(cand,rep_lookup[cand.rep_name],wavelengths,Xg[tr],y[tr],stable_seed(seed,cand.key,j))
            diag=_extract_ground_elasticnet_diagnostic(fit)
            if diag is not None:
                d=dict(diag); d["inner_fold"]=int(j); d["candidate_key"]=cand.key; d["family"]=cand.family; elastic_audit.append(d)
            pred=predict_ground_candidate(fit,Xg[va])
            if not np.all(np.isfinite(pred)):
                return float("inf"),float("inf"),[],elastic_audit,"invalid","NonFinitePrediction",f"non-finite prediction in inner fold {j}"
            rmses.append(float(np.sqrt(np.mean((pred-y[va])**2))))
        except ElasticNetConvergenceError as exc:
            d=dict(exc.diagnostic); d["inner_fold"]=int(j); d["candidate_key"]=cand.key; d["family"]=cand.family; elastic_audit.append(d)
            return float("inf"),float("inf"),[],elastic_audit,"invalid","ElasticNetConvergenceError",str(exc)
        except (np.linalg.LinAlgError, FloatingPointError) as exc:
            return float("inf"),float("inf"),[],elastic_audit,"invalid",type(exc).__name__,str(exc)
    if len(rmses)!=len(inner):
        return float("inf"),float("inf"),[],elastic_audit,"invalid","IncompleteInnerCV","not all inner folds produced predictions"
    return float(np.mean(rmses)),float(np.std(rmses,ddof=1)/math.sqrt(len(rmses))) if len(rmses)>1 else 0.0,rmses,elastic_audit,"valid","",""

def screen_representations_branch(reps: List[Representation], wavelengths: np.ndarray,
                                  Xw: np.ndarray, Xg: np.ndarray, y: np.ndarray,
                                  inner: List[Tuple[np.ndarray,np.ndarray]], seed: int,
                                  top_k: int, smoke: bool, branch: str) -> Tuple[List[Representation],pd.DataFrame]:







    branch=str(branch).lower(); rep_lookup={r.name:r for r in reps}; rows=[]
    if branch not in {"whole","transfer","ground"}: raise ValueError(branch)
    for rep in reps:
        cands=_screen_proxy_candidates(rep,smoke,branch)
        candidate_scores=[]
        for cand in cands:
            if branch == "ground":
                mean,se,_,_,status,failure_type,failure_message=cv_ground_candidate(cand,rep_lookup,wavelengths,Xg,y,inner,stable_seed(seed,branch,"rep_screen",cand.family))
            else:
                mean,se,_,_,status,failure_type,failure_message=cv_candidate(cand,rep_lookup,wavelengths,Xw,Xg,y,inner,stable_seed(seed,branch,"rep_screen",cand.family))
            candidate_scores.append((float(mean),cand.family,cand.params_json,status,failure_type,failure_message))
        all_valid=(len(candidate_scores)==4 and all(np.isfinite(x[0]) and x[3]=="valid" for x in candidate_scores))
        if all_valid:
            vals=np.asarray([x[0] for x in candidate_scores],float)
            screen_score=float(np.median(vals))
            best=min(candidate_scores,key=lambda x:(x[0],x[1],x[2]))
            best_rmse,best_family,best_params,_,_,_=best
        else:
            screen_score=float("inf"); best_rmse=float("inf"); best_family="none"; best_params="{}"
        proxy_json=json.dumps([
            {"family":fam,"params_json":params,"mean_inner_rmse":rmse,"fit_status":status,"failure_type":failure_type,"failure_message":failure_message}
            for rmse,fam,params,status,failure_type,failure_message in candidate_scores
        ],separators=(",",":"))
        rows.append({"branch":branch,"representation":rep.name,"screen_rmse":screen_score,
                     "screen_score_definition":"median_mean_inner_rmse_across_4_fixed_family_spanning_proxies",
                     "screen_best_proxy_rmse":best_rmse,"screen_best_family":best_family,
                     "screen_best_params_json":best_params,"screen_proxy_scores_json":proxy_json,
                     "screen_candidate_count":len(candidate_scores),"all_screen_proxies_valid":bool(all_valid),
                     "min_nm":rep.min_nm,"max_nm":rep.max_nm,"preprocess":rep.preprocess,
                     "sg_window":rep.sg_window,"sg_poly":rep.sg_poly})
    df=pd.DataFrame(rows).sort_values(["screen_rmse","representation"]).reset_index(drop=True)
    valid=df[np.isfinite(df["screen_rmse"])].copy()
    if len(valid)<top_k:
        raise RuntimeError(f"Branch {branch} has only {len(valid)} valid representations; {top_k} required")
    selected=[rep_lookup[n] for n in valid.head(top_k)["representation"]]
    return selected,df

def tune_ground_candidates(cands: List[Candidate], reps: List[Representation], wavelengths: np.ndarray,
                           Xg: np.ndarray, y: np.ndarray, inner: List[Tuple[np.ndarray,np.ndarray]],
                           seed: int) -> Tuple[Candidate,pd.DataFrame,Dict[str,Candidate]]:
    lookup={r.name:r for r in reps}; rows=[]
    for cand in cands:
        mean,se,folds,elastic_audit,fit_status,failure_type,failure_message=cv_ground_candidate(cand,lookup,wavelengths,Xg,y,inner,seed)
        ea=elastic_audit
        rows.append({"candidate_key":cand.key,"family":cand.family,"representation":cand.rep_name,"params_json":cand.params_json,
                     "complexity":cand.complexity,"mean_rmse":mean,"se_rmse":se,"fold_rmse_json":json.dumps(folds),
                     "fit_status":fit_status,"failure_type":failure_type,"failure_message":failure_message,
                     "elasticnet_candidate":bool(cand.family=="ground_elastic"),
                     "elasticnet_all_inner_folds_converged":(bool(len(ea)==len(inner) and all(bool(d.get("converged",False)) for d in ea)) if cand.family=="ground_elastic" else np.nan),
                     "elasticnet_retry_count_total":(int(sum(max(0,int(d.get("attempts_used",1))-1) for d in ea)) if ea else 0),
                     "elasticnet_max_n_iter":(int(max([int(d.get("final_n_iter",0)) for d in ea],default=0)) if ea else 0),
                     "elasticnet_max_dual_gap":(float(max([float(d.get("final_dual_gap",0.0)) for d in ea],default=0.0)) if ea else np.nan),
                     "elasticnet_audit_json":json.dumps(ea,separators=(",",":"))})
    df=pd.DataFrame(rows); chosen_row=select_min_inner_rmse(df); cmap={c.key:c for c in cands}; chosen=cmap[str(chosen_row["candidate_key"])]
    family_best={}
    for family,g in df.groupby("family"):
        if not np.isfinite(g["mean_rmse"]).any(): continue
        rr=select_min_inner_rmse(g); family_best[family]=cmap[str(rr["candidate_key"])]
    return chosen,df,family_best

def tune_candidates(cands: List[Candidate], reps: List[Representation], wavelengths: np.ndarray, Xw: np.ndarray, Xg: np.ndarray, y: np.ndarray,
                    inner: List[Tuple[np.ndarray,np.ndarray]], seed: int) -> Tuple[Candidate,pd.DataFrame,Dict[str,Candidate]]:
    lookup={r.name:r for r in reps}; rows=[]
    for cand in cands:
        mean,se,folds,elastic_audit,fit_status,failure_type,failure_message=cv_candidate(cand,lookup,wavelengths,Xw,Xg,y,inner,seed)
        ea = elastic_audit
        rows.append({"candidate_key":cand.key,"family":cand.family,"representation":cand.rep_name,"params_json":cand.params_json,
                     "complexity":cand.complexity,"mean_rmse":mean,"se_rmse":se,"fold_rmse_json":json.dumps(folds),
                     "fit_status":fit_status,"failure_type":failure_type,"failure_message":failure_message,
                     "elasticnet_candidate":bool(cand.family in {"whole_elastic","epo_elastic","pooled_elastic","ground_elastic"}),
                     "elasticnet_all_inner_folds_converged":(bool(len(ea)==len(inner) and all(bool(d.get("converged",False)) for d in ea)) if cand.family in {"whole_elastic","epo_elastic","pooled_elastic","ground_elastic"} else np.nan),
                     "elasticnet_retry_count_total":(int(sum(max(0,int(d.get("attempts_used",1))-1) for d in ea)) if ea else 0),
                     "elasticnet_max_n_iter":(int(max([int(d.get("final_n_iter",0)) for d in ea],default=0)) if ea else 0),
                     "elasticnet_max_dual_gap":(float(max([float(d.get("final_dual_gap",0.0)) for d in ea],default=0.0)) if ea else np.nan),
                     "elasticnet_audit_json":json.dumps(ea,separators=(",",":"))})
    df=pd.DataFrame(rows)
    chosen_row=select_min_inner_rmse(df)
    cmap={c.key:c for c in cands}
    chosen=cmap[str(chosen_row["candidate_key"])]
    family_best={}
    for family,g in df.groupby("family"):
        if not np.isfinite(g["mean_rmse"]).any():


            continue
        rr=select_min_inner_rmse(g); family_best[family]=cmap[str(rr["candidate_key"])]
    return chosen,df,family_best








def summarize_oof(oof:pd.DataFrame)->pd.DataFrame:
    rows=[]
    for (target,method,repeat),g in oof.groupby(["target","method","outer_repeat"]):
        m=regression_metrics(g["y_true"].to_numpy(float),g["y_pred"].to_numpy(float))
        rows.append({"target":target,"method":method,"outer_repeat":repeat,**m})
    rep=pd.DataFrame(rows)
    out=[]
    for (target,method),g in rep.groupby(["target","method"]):
        row={"target":target,"method":method,"n_outer_repeats":len(g)}
        for metric in ["r2","rmse","mae","bias","rpd"]:
            row[f"median_{metric}"]=float(g[metric].median()); row[f"q25_{metric}"]=float(g[metric].quantile(.25)); row[f"q75_{metric}"]=float(g[metric].quantile(.75))
        out.append(row)
    return pd.DataFrame(out),rep


def _inference_group_column(df: pd.DataFrame) -> str:
    return "inference_group" if "inference_group" in df.columns else "sample_id"

def paired_cluster_bootstrap(oof:pd.DataFrame,target:str,a:str,b:str,n_boot:int,seed:int)->Dict[str,float]:
    df=oof[(oof.target==target)&(oof.method.isin([a,b]))]
    cluster_col=_inference_group_column(df)
    clusters=np.array(sorted(df[cluster_col].astype(str).unique()),dtype=str)
    per_cluster=[]
    for cluster in clusters:
        ga=df[df[cluster_col].astype(str)==cluster]
        aa=ga[ga.method==a]; bb=ga[ga.method==b]
        ma=aa.merge(bb,on=["sample_id","outer_repeat","outer_fold",cluster_col] if cluster_col!="sample_id" else ["sample_id","outer_repeat","outer_fold"],suffixes=("_a","_b"))
        if len(ma):
            per_cluster.append(float(np.mean((ma.y_pred_a-ma.y_true_a)**2-(ma.y_pred_b-ma.y_true_b)**2)))
    vals=np.asarray(per_cluster,float)
    if len(vals)==0: raise RuntimeError(f"No paired inference clusters for {target}: {a} vs {b}")
    observed=float(vals.mean())
    rng=np.random.default_rng(seed); boots=np.empty(n_boot,float)
    for i in range(n_boot): boots[i]=float(rng.choice(vals,size=len(vals),replace=True).mean())
    return {"target":target,"method_a":a,"method_b":b,"observed_mse_improvement_b_over_a":observed,"mse_improvement_b_over_a_mean":observed,"ci_low":float(np.quantile(boots,.025)),"ci_high":float(np.quantile(boots,.975)),"bootstrap_repeats":int(n_boot),"inference_cluster_column":cluster_col,"n_inference_clusters":int(len(vals))}

def sample_cluster_signflip(oof:pd.DataFrame,target:str,a:str,b:str,n_perm:int,seed:int)->Dict[str,float]:
    df=oof[(oof.target==target)&(oof.method.isin([a,b]))]
    cluster_col=_inference_group_column(df)
    clusters=np.array(sorted(df[cluster_col].astype(str).unique()),dtype=str); d=[]
    for cluster in clusters:
        ga=df[df[cluster_col].astype(str)==cluster]
        aa=ga[ga.method==a]; bb=ga[ga.method==b]
        merge_keys=["sample_id","outer_repeat","outer_fold"]
        if cluster_col!="sample_id": merge_keys.append(cluster_col)
        ma=aa.merge(bb,on=merge_keys,suffixes=("_a","_b"))
        d.append(float(np.mean((ma.y_pred_a-ma.y_true_a)**2-(ma.y_pred_b-ma.y_true_b)**2)))
    d=np.asarray(d,float); obs=float(np.mean(d)); rng=np.random.default_rng(seed); null=np.empty(n_perm)
    for i in range(n_perm): null[i]=float(np.mean(d*rng.choice([-1.0,1.0],size=len(d))))
    p=float((1+np.sum(null>=obs))/(n_perm+1))
    return {"target":target,"method_a":a,"method_b":b,"observed_mean_mse_improvement_b_over_a":obs,"one_sided_p":p,"n_inference_clusters":len(d),"inference_cluster_column":cluster_col,"permutations":n_perm}

def holm(pvals:Sequence[float])->np.ndarray:
    p=np.asarray(pvals,float); order=np.argsort(p); out=np.empty_like(p); running=0.0; m=len(p)
    for rank,idx in enumerate(order):
        val=min(1.0,(m-rank)*p[idx]); running=max(running,val); out[idx]=running
    return out






def analyze_target(target:str,data:PairedSpectralData,args:argparse.Namespace,out_dir:Path,cv_groups:Optional[np.ndarray])->Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    y=data.chemistry.set_index("sample_id").loc[data.sample_ids,target].to_numpy(float)
    Xw=data.whole_mean; Xg=data.ground_mean; waves=data.wavelengths
    inference_groups=np.asarray(cv_groups,dtype=str) if cv_groups is not None else data.sample_ids.astype(str)
    reps=representation_space(args.profile); rep_lookup={r.name:r for r in reps}
    outer=make_regression_stratified_splits(y,args.outer_splits,args.outer_repeats,stable_seed(args.seed,target),groups=cv_groups)
    oof_rows=[]; tune_rows=[]; screen_rows=[]; selection_rows=[]; numerical_rows=[]; replicate_prediction_rows=[]
    smoke=args.profile=="smoke"
    for repeat,fold,tr,te in outer:
        fold_seed=stable_seed(args.seed,target,repeat,fold)
        inner_groups=cv_groups[tr] if cv_groups is not None else None
        inner=[(a,b) for _,_,a,b in make_regression_stratified_splits(y[tr],args.inner_splits,1,fold_seed,groups=inner_groups)]


        whole_reps,whole_screen=screen_representations_branch(reps,waves,Xw[tr],Xg[tr],y[tr],inner,fold_seed,args.top_representations_per_branch,smoke,"whole")
        transfer_reps,transfer_screen=screen_representations_branch(reps,waves,Xw[tr],Xg[tr],y[tr],inner,fold_seed,args.top_representations_per_branch,smoke,"transfer")
        ground_reps,ground_screen=screen_representations_branch(reps,waves,Xw[tr],Xg[tr],y[tr],inner,fold_seed,args.top_representations_per_branch,smoke,"ground")
        for screen in [whole_screen,transfer_screen,ground_screen]:
            screen["target"]=target; screen["outer_repeat"]=repeat; screen["outer_fold"]=fold; screen_rows.append(screen)
        whole_cands=[]; transfer_cands=[]; ground_cands=[]
        for rep in whole_reps: whole_cands.extend(_candidate_grids(rep,args.profile,transfer=False))
        for rep in transfer_reps: transfer_cands.extend(_candidate_grids(rep,args.profile,transfer=True))
        for rep in ground_reps: ground_cands.extend(_ground_candidate_grids(rep,args.profile))
        whole_best,whole_trace,whole_family=tune_candidates(whole_cands,whole_reps,waves,Xw[tr],Xg[tr],y[tr],inner,fold_seed)
        trans_best,trans_trace,trans_family=tune_candidates(transfer_cands,transfer_reps,waves,Xw[tr],Xg[tr],y[tr],inner,fold_seed)
        ground_best,ground_trace,ground_family=tune_ground_candidates(ground_cands,ground_reps,waves,Xg[tr],y[tr],inner,fold_seed)
        for typ,trace in [("whole",whole_trace),("transfer",trans_trace),("ground_reference",ground_trace)]:
            trace=trace.copy(); trace["candidate_group"]=typ; trace["target"]=target; trace["outer_repeat"]=repeat; trace["outer_fold"]=fold; tune_rows.append(trace)

        fitted:Dict[str,Tuple[Candidate,Dict[str,Any]]]={}
        fitted["nested_strong_whole"]=(whole_best,fit_candidate(whole_best,rep_lookup[whole_best.rep_name],waves,Xw[tr],Xg[tr],y[tr],fold_seed))
        fitted["nested_ground_assisted"]=(trans_best,fit_candidate(trans_best,rep_lookup[trans_best.rep_name],waves,Xw[tr],Xg[tr],y[tr],fold_seed))
        for fam,cand in whole_family.items(): fitted[fam]=(cand,fit_candidate(cand,rep_lookup[cand.rep_name],waves,Xw[tr],Xg[tr],y[tr],fold_seed))
        for fam,cand in trans_family.items(): fitted[fam]=(cand,fit_candidate(cand,rep_lookup[cand.rep_name],waves,Xw[tr],Xg[tr],y[tr],fold_seed))
        for method,(cand,fitobj) in fitted.items():
            diag=_extract_elasticnet_diagnostic(fitobj)
            if diag is not None:
                row=dict(diag); row.update({"target":target,"outer_repeat":int(repeat),"outer_fold":int(fold),
                                            "method":method,"candidate_key":cand.key,"family":cand.family,
                                            "representation":cand.rep_name,"params_json":cand.params_json})
                numerical_rows.append(row)
        ground=fit_ground_candidate(ground_best,rep_lookup[ground_best.rep_name],waves,Xg[tr],y[tr],fold_seed)
        gdiag=_extract_ground_elasticnet_diagnostic(ground)
        if gdiag is not None:
            row=dict(gdiag); row.update({"target":target,"outer_repeat":int(repeat),"outer_fold":int(fold),
                                        "method":"ground_reference","candidate_key":ground_best.key,"family":ground_best.family,
                                        "representation":ground_best.rep_name,"params_json":ground_best.params_json})
            numerical_rows.append(row)

        for method,(cand,fit) in fitted.items():
            pred=predict_candidate(fit,Xw[te])
            for j,idx in enumerate(te):
                oof_rows.append({"target":target,"method":method,"sample_id":data.sample_ids[idx],"inference_group":inference_groups[idx],"outer_repeat":repeat,"outer_fold":fold,"y_true":y[idx],"y_pred":pred[j],"selected_family":cand.family,"selected_representation":cand.rep_name,"selected_params_json":cand.params_json})
        gp=predict_ground_candidate(ground,Xg[te])
        for j,idx in enumerate(te):
            oof_rows.append({"target":target,"method":"ground_reference","sample_id":data.sample_ids[idx],"inference_group":inference_groups[idx],"outer_repeat":repeat,"outer_fold":fold,"y_true":y[idx],"y_pred":gp[j],"selected_family":ground_best.family,"selected_representation":ground_best.rep_name,"selected_params_json":ground_best.params_json})

        inner_manifest=[]
        for inner_fold,(itr,iva) in enumerate(inner):
            inner_manifest.append({
                "inner_fold":int(inner_fold),
                "train_sample_ids":data.sample_ids[tr][itr].tolist(),
                "validation_sample_ids":data.sample_ids[tr][iva].tolist(),
            })
        selection_rows.append({"target":target,"outer_repeat":repeat,"outer_fold":fold,
                               "whole_selected_candidate_key":whole_best.key,"whole_selected_family":whole_best.family,"whole_selected_representation":whole_best.rep_name,"whole_selected_params_json":whole_best.params_json,
                               "transfer_selected_candidate_key":trans_best.key,"transfer_selected_family":trans_best.family,"transfer_selected_representation":trans_best.rep_name,"transfer_selected_params_json":trans_best.params_json,
                               "ground_selected_candidate_key":ground_best.key,"ground_selected_family":ground_best.family,"ground_selected_representation":ground_best.rep_name,"ground_selected_params_json":ground_best.params_json,
                               "whole_screen_selected_representations_json":json.dumps([r.name for r in whole_reps]),
                               "transfer_screen_selected_representations_json":json.dumps([r.name for r in transfer_reps]),
                               "ground_screen_selected_representations_json":json.dumps([r.name for r in ground_reps]),
                               "representation_candidates_per_branch":len(reps),
                               "representations_retained_per_branch":args.top_representations_per_branch,
                               "outer_train_sample_ids_json":json.dumps(data.sample_ids[tr].tolist()),
                               "outer_test_sample_ids_json":json.dumps(data.sample_ids[te].tolist()),
                               "inner_split_sample_ids_json":json.dumps(inner_manifest)})
        print(f"[{target}] repeat {repeat+1}/{args.outer_repeats} fold {fold+1}/{args.outer_splits}: whole={whole_best.family} transfer={trans_best.family} ground={ground_best.family}")
    return pd.DataFrame(oof_rows),pd.concat(tune_rows,ignore_index=True),pd.concat(screen_rows,ignore_index=True),pd.DataFrame(selection_rows),pd.DataFrame(numerical_rows),pd.DataFrame(replicate_prediction_rows)



def windows_path_preflight(package_root:Path, results_dir:Optional[Path]=None, margin:int=240)->Dict[str,Any]:







    package_root=Path(package_root).resolve()
    results_dir=Path(results_dir).resolve() if results_dir is not None else package_root/'v13_2_results'
    representatives=[
        package_root/'experimental_code'/'sorghum_v13_1_frozen.py',
        package_root/'experimental_data'/'whole_grain_replicate_spectra.csv.gz',
        results_dir/'00_audit'/'elasticnet_inner_cv_convergence_audit.csv.gz',
        results_dir/'04_inference'/'ground_assisted_vs_strong_whole_cluster_bootstrap.csv',
        results_dir/'02_selection'/'outer_fold_selected_workflows.csv',
    ]
    longest=max(representatives,key=lambda x:len(str(x)))
    report={"platform_is_windows":os.name=='nt',"margin_characters":int(margin),
            "longest_representative_path":str(longest),"longest_representative_length":len(str(longest)),
            "safe":True}
    if os.name=='nt' and len(str(longest))>=margin:
        report['safe']=False
        raise RuntimeError(
            "WINDOWS_PATH_UNSAFE: the current package location is too deep for reliable formal V13.2.0-98 analysis "
            f"(representative path length={len(str(longest))}, safety margin={margin}). "
            "Move/extract this package to a short directory such as C:\\NIR_V13_1 and run again. "
            "Do not place the ZIP inside another same-named folder."
        )
    return report

def validate_reference_metadata(data_dir:Path)->Dict[str,Any]:
    path=data_dir/"reference_method_metadata_template.csv"
    report={"file_present":path.exists(),"complete":False,"missing":[],"structure_valid":False}
    if not path.exists(): return report
    df=pd.read_csv(path).fillna("")
    if "target" not in df.columns: raise ValueError("reference_method_metadata_template.csv must contain target")
    if df["target"].astype(str).duplicated().any(): raise ValueError("Duplicate target rows in reference method metadata")
    if set(df["target"].astype(str))!=set(TARGETS): raise ValueError(f"Reference method metadata targets must be exactly {TARGETS}")
    req=["standard_method_or_protocol","analytical_replicates","reporting_basis","reference_instrument_or_apparatus","laboratory_or_institution"]
    missing_cols=[c for c in req if c not in df.columns]
    if missing_cols: raise ValueError(f"Reference method metadata missing required columns: {missing_cols}")
    for _,r in df.iterrows():
        target=str(r["target"])
        for c in req:
            if str(r.get(c,"")).strip()=="": report["missing"].append(f"{target}:{c}")
        rep=str(r.get("analytical_replicates","")).strip()
        if rep:
            try:
                if int(float(rep))<1 or abs(float(rep)-int(float(rep)))>1e-12: raise ValueError
            except ValueError:
                raise ValueError(f"analytical_replicates must be a positive integer for {target}")
    report["structure_valid"]=True
    report["complete"]=len(report["missing"])==0
    return report


def environment_manifest(args,data_dir)->Dict[str,Any]:
    return {"framework_version":FRAMEWORK_VERSION,"script_version":SCRIPT_VERSION,"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"scikit_learn":sklearn.__version__,
            "seed":args.seed,"seed_is_selectable":False,"profile":args.profile,"outer_splits":args.outer_splits,"outer_repeats":args.outer_repeats,"inner_splits":args.inner_splits,
            "development_status":"development cohort; independent external validation not performed by this script",
            "fairness_correction":{
                "common_representation_space_for_all_branches":True,
                "branch_specific_representation_screening":True,"family_spanning_proxy_screening":True,"screen_proxy_count_per_branch":4,"screen_score_aggregation":"median",
                "same_inner_splits_for_all_branches":True,"same_proxy_count_and_score_aggregation_for_all_branches":True,
                "same_representations_retained_per_branch":int(args.top_representations_per_branch),
                "ground_preprocessing_fit_exclusively_on_ground_training":True,
                "ground_base_model_families_mirror_whole_base_model_families":True,
                "outer_test_used_for_selection":False
            },
            "input_sha256":{f:sha256_file(data_dir/f) for f in ["whole_grain_replicate_spectra.csv.gz","ground_grain_replicate_spectra.csv.gz","reference_chemistry.csv","sample_metadata_template.csv","reference_method_metadata_template.csv"] if (data_dir/f).exists()}}


def build_parser()->argparse.ArgumentParser:
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir",type=Path,default=Path(__file__).resolve().parents[1]/"experimental_data")
    p.add_argument("--results-dir",type=Path,default=Path(__file__).resolve().parents[1]/"v13_2_results")
    p.add_argument("--targets",nargs="+",default=TARGETS,choices=TARGETS)
    p.add_argument("--profile",choices=["smoke","publication","extended"],default="publication")
    p.add_argument("--seed",type=int,default=DEFAULT_SEED,help="Fixed reproducibility seed; not selected by performance.")
    p.add_argument("--outer-splits",type=int,default=5); p.add_argument("--outer-repeats",type=int,default=10); p.add_argument("--inner-splits",type=int,default=4)
    p.add_argument("--top-representations-per-branch",type=int,default=2,help="Same number of independently screened representations retained by whole, transfer, and ground branches.")
    p.add_argument("--bootstrap-repeats",type=int,default=2000); p.add_argument("--signflip-repeats",type=int,default=20000)
    p.add_argument("--grid-step-nm",type=float,default=2.0)
    p.add_argument("--smoke-sample-limit",type=int,default=0,help="Smoke-only integration-test sample cap; ignored in publication profile.")
    return p


def enforce_publication_design_lock(args: argparse.Namespace) -> None:

    if args.profile != "publication":
        return
    expected = {
        "seed": DEFAULT_SEED,
        "outer_splits": 5,
        "outer_repeats": 10,
        "inner_splits": 4,
        "top_representations_per_branch": 2,
        "bootstrap_repeats": 2000,
        "signflip_repeats": 20000,
        "grid_step_nm": 2.0,
    }
    bad = []
    for key, value in expected.items():
        actual = getattr(args, key)
        if isinstance(value, float):
            ok = abs(float(actual)-value) < 1e-12
        else:
            ok = actual == value
        if not ok:
            bad.append(f"{key}={actual!r} (locked={value!r})")
    if bad:
        raise RuntimeError(
            "PUBLICATION_DESIGN_LOCK_VIOLATION: formal V13.2.0-98 settings are prespecified and cannot "
            "be changed in publication profile. Use smoke/extended only for software testing or "
            "clearly labeled exploratory work. Violations: " + "; ".join(bad)
        )


def main()->None:
    args=build_parser().parse_args(); enforce_publication_design_lock(args); windows_path_preflight(Path(__file__).resolve().parents[1],args.results_dir); out=ensure_dir(args.results_dir); data=load_data(args.data_dir,args.grid_step_nm,400.0,2498.0)
    cv_groups,metadata_report=load_metadata_groups(args.data_dir,data.sample_ids)
    if args.profile=="smoke" and int(args.smoke_sample_limit)>0:
        n=min(int(args.smoke_sample_limit),len(data.sample_ids))
        if n < max(8,args.outer_splits*args.inner_splits):
            raise ValueError("smoke_sample_limit is too small for the requested nested CV")
        ids=data.sample_ids[:n].copy()
        chem=data.chemistry.set_index("sample_id").loc[ids].reset_index()
        data=PairedSpectralData(ids,data.wavelengths,data.whole_replicates[:n].copy(),data.ground_replicates[:n].copy(),chem)
        if cv_groups is not None: cv_groups=cv_groups[:n].copy()
        metadata_report=dict(metadata_report); metadata_report["smoke_sample_limit_applied"]=n
    ensure_dir(out/"00_audit"); save_json(environment_manifest(args,args.data_dir),out/"00_audit"/"run_manifest.json"); save_json(metadata_report,out/"00_audit"/"sample_provenance_audit.json"); save_json(validate_reference_metadata(args.data_dir),out/"00_audit"/"reference_method_metadata_audit.json")
    reps_for_audit=representation_space(args.profile)
    fairness_audit={
        "status":"passed",
        "framework_version":FRAMEWORK_VERSION,"core_script_version":SCRIPT_VERSION,
        "profile":args.profile,
        "common_representation_space_count":len(reps_for_audit),
        "whole_representation_space_count":len(reps_for_audit),
        "transfer_representation_space_count":len(reps_for_audit),
        "ground_representation_space_count":len(reps_for_audit),
        "representations_retained_per_branch":int(args.top_representations_per_branch),
        "representation_screening_rule":"Each branch evaluates every representation with exactly four fixed family-spanning proxies and ranks by the median proxy mean inner-CV RMSE.",
        "screen_proxy_count_per_branch":4,
        "screen_score_aggregation":"median of four proxy mean inner-CV RMSE values",
        "whole_representation_screen_families":["whole_pls","whole_elastic","whole_svr","whole_krr"],
        "transfer_representation_screen_families":["epo_pls","epo_elastic","pooled_elastic","dipls"],
        "ground_representation_screen_families":["ground_pls","ground_elastic","ground_svr","ground_krr"],
        "ground_full_model_families":["ground_pls","ground_ridge","ground_elastic","ground_ipls","ground_svr","ground_krr","ground_pls_krr","ground_extratrees"],
        "whole_full_model_families":["whole_pls","whole_ridge","whole_elastic","whole_ipls","whole_svr","whole_krr","whole_pls_krr","whole_extratrees"],
        "transfer_full_model_families":["epo_pls","epo_elastic","pooled_elastic","dipls"],
        "same_outer_and_inner_splits":True,
        "outer_test_used_for_screening_or_tuning":False,
        "ground_branch_depends_on_whole_selected_representation":False,
        "transfer_branch_depends_on_whole_selected_representation":False
    }
    save_json(fairness_audit,out/"00_audit"/"model_space_fairness_audit.json")
    all_oof=[];all_tune=[];all_screen=[];all_sel=[];all_num=[];all_rep_pred=[]
    for target in args.targets:
        print(f"[{target}] starting V13.2.0-98 frozen fully nested boundary analysis")
        o,t,scr,sel,num,rp=analyze_target(target,data,args,ensure_dir(out/target),cv_groups); all_oof.append(o);all_tune.append(t);all_screen.append(scr);all_sel.append(sel);all_num.append(num);all_rep_pred.append(rp)
    oof=pd.concat(all_oof,ignore_index=True); tune=pd.concat(all_tune,ignore_index=True); screen=pd.concat(all_screen,ignore_index=True); sel=pd.concat(all_sel,ignore_index=True); numerical=pd.concat(all_num,ignore_index=True) if any(len(x) for x in all_num) else pd.DataFrame(); replicate_predictions=pd.concat(all_rep_pred,ignore_index=True) if any(len(x) for x in all_rep_pred) else pd.DataFrame()
    ensure_dir(out/"01_predictions"); oof.to_csv(out/"01_predictions"/"all_targets_outer_oof_predictions.csv",index=False)
    ensure_dir(out/"06_robustness")
    if not replicate_predictions.empty: replicate_predictions.to_csv(out/"06_robustness"/"technical_replicate_oof_predictions.csv",index=False)
    ensure_dir(out/"02_selection"); tune.to_csv(out/"02_selection"/"inner_cv_candidate_trace.csv.gz",index=False,compression="gzip"); screen.to_csv(out/"02_selection"/"representation_screening_trace.csv",index=False); sel.to_csv(out/"02_selection"/"outer_fold_selected_workflows.csv",index=False)
    ensure_dir(out/"00_audit")
    elastic_inner=tune[tune.get("elasticnet_candidate",False)==True].copy() if "elasticnet_candidate" in tune.columns else pd.DataFrame()
    if not elastic_inner.empty:
        elastic_inner.to_csv(out/"00_audit"/"elasticnet_inner_cv_convergence_audit.csv.gz",index=False,compression="gzip")
    if not numerical.empty:
        numerical.to_csv(out/"00_audit"/"elasticnet_selected_and_family_fit_audit.csv",index=False)
    numerical_summary={
        "policy":"ElasticNet warnings are captured; deterministic retry is attempted; unresolved candidates are invalid and cannot enter selection.",
        "inner_elastic_candidates_evaluated":int(len(elastic_inner)),
        "inner_elastic_candidates_all_folds_converged":int(elastic_inner["elasticnet_all_inner_folds_converged"].astype("boolean").fillna(False).sum()) if not elastic_inner.empty else 0,
        "inner_elastic_candidates_invalid_or_nonconverged":int((~elastic_inner["elasticnet_all_inner_folds_converged"].astype("boolean").fillna(False)).sum()) if not elastic_inner.empty else 0,
        "final_elastic_fits_recorded":int(len(numerical)),
        "final_elastic_fits_converged":int(numerical["converged"].fillna(False).sum()) if not numerical.empty else 0,
        "unresolved_selected_final_fits":int((~numerical["converged"].fillna(False)).sum()) if not numerical.empty else 0,
    }
    save_json(numerical_summary,out/"00_audit"/"numerical_convergence_summary.json")
    if numerical_summary["unresolved_selected_final_fits"] != 0:
        raise RuntimeError("Strict numerical QA failed: at least one selected/family ElasticNet final fit did not converge")
    split_manifest=oof[oof.method=="nested_strong_whole"][["target","sample_id","inference_group","outer_repeat","outer_fold"]].drop_duplicates().sort_values(["target","outer_repeat","outer_fold","sample_id"])
    split_manifest.to_csv(out/"02_selection"/"outer_test_split_manifest.csv",index=False)
    perf,repeat_perf=summarize_oof(oof); ensure_dir(out/"03_performance"); perf.to_csv(out/"03_performance"/"method_performance_summary.csv",index=False); repeat_perf.to_csv(out/"03_performance"/"outer_repeat_performance.csv",index=False)

    inf=[]; boots=[]
    for target in args.targets:
        inf.append(sample_cluster_signflip(oof,target,"nested_strong_whole","nested_ground_assisted",args.signflip_repeats,stable_seed(args.seed,"signflip",target)))
        boots.append(paired_cluster_bootstrap(oof,target,"nested_strong_whole","nested_ground_assisted",args.bootstrap_repeats,stable_seed(args.seed,"bootstrap",target)))
    infdf=pd.DataFrame(inf); infdf["role"]=infdf.target.map(lambda t:"primary" if t==PRIMARY_TARGET else ("secondary" if t==SECONDARY_TARGET else "supportive"))

    mask=infdf.target!=PRIMARY_TARGET
    infdf["holm_adjusted_p_nonprimary"]=np.nan
    if mask.sum(): infdf.loc[mask,"holm_adjusted_p_nonprimary"]=holm(infdf.loc[mask,"one_sided_p"].to_numpy())
    ensure_dir(out/"04_inference"); infdf.to_csv(out/"04_inference"/"ground_assisted_vs_strong_whole_signflip.csv",index=False); pd.DataFrame(boots).to_csv(out/"04_inference"/"ground_assisted_vs_strong_whole_cluster_bootstrap.csv",index=False)

    freq=[]
    for col,label in [("whole_selected_family","whole"),("transfer_selected_family","transfer"),("ground_selected_family","ground_reference")]:
        g=sel.groupby(["target",col]).size().reset_index(name="n_outer_folds"); g["selection_group"]=label; g=g.rename(columns={col:"family"}); freq.append(g)
    pd.concat(freq,ignore_index=True).to_csv(out/"03_performance"/"nested_method_selection_frequency.csv",index=False)
    summary={"script_version":SCRIPT_VERSION,"coded_samples":len(data.sample_ids),"targets":args.targets,"primary_target":PRIMARY_TARGET,
             "primary_comparison":"nested_ground_assisted vs nested_strong_whole, selected wholly inside each outer-training set",
             "seed_selection":False,"outer_test_used_for_model_selection":False,"external_validation_performed":False,"strict_numerical_convergence_policy":True,
             "important_limitation":"The method space was informed by previous analyses of this development cohort. Independent external validation is required for confirmatory deployment claims."}
    save_json(summary,out/"ANALYSIS_SUMMARY.json")
    print(f"V13.2.0-98 core analysis complete. Results: {out}")

if __name__=="__main__": main()
