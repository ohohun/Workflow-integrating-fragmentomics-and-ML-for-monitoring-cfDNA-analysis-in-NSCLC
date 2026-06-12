from __future__ import annotations

import hashlib
import hmac
import io
import os
import re
import json
import zipfile
import secrets
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Optional, Any

import pandas as pd
import numpy as np
from fastapi import Depends, FastAPI, HTTPException, Response, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.security import APIKeyCookie

# --- Configuration ---
BUCKET = "cfdna-hospital-runtime-36ccdf91-f5b0-4219-87e"
USE_LOCAL_DATA = os.getenv("USE_LOCAL_DATA", "1") == "1"

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

# Standard paths (Relative to DATA_DIR)
MASTER_CSV = "master/master_hospital.csv"
ML_REGISTRY_CSV = "ml/ml_registry_v4.csv"
ML_TRAJECTORIES_CSV = "ml/ml_trajectories.csv"
ENDMOTIF_CSV = "features/endmotif_merged_sorted.csv"
FRAG_BINS_CSV = "features/fraglen_merged_sorted.csv"
REPORTS_DIR = "reports"

_cache = {}

app = FastAPI(title="cfDNA Monitor API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

cookie_scheme = APIKeyCookie(name="session_id", auto_error=False)

# --- Auth + RBAC + Admin user management ---
SESSION_TTL_HOURS = int(os.getenv("SESSION_TTL_HOURS", "8"))
ROLE_RANK = {"viewer": 0, "analyst": 1, "admin": 2}
USER_STORE_PATH = DATA_DIR / "users.json"
AUDIT_LOG_PATH = DATA_DIR / "audit_log.json"

def _hash_password(password: str, salt_hex: Optional[str] = None):
    salt = bytes.fromhex(salt_hex) if salt_hex else secrets.token_bytes(32)
    dk = hashlib.pbkdf2_hmac("sha256", password.encode(), salt, 260_000)
    return salt.hex(), dk.hex()

def _verify_password(password: str, salt_hex: str, dk_hex: str) -> bool:
    _, new_dk = _hash_password(password, salt_hex)
    return hmac.compare_digest(new_dk, dk_hex)

def _make_user(password: str, role: str, display: str) -> dict:
    salt, dk = _hash_password(password)
    now = datetime.utcnow().isoformat()
    return {
        "salt": salt,
        "dk": dk,
        "role": role,
        "display": display,
        "is_active": True,
        "created_at": now,
        "updated_at": now,
        "last_login": None,
    }

DEFAULT_USERS: Dict[str, dict] = {
    "admin": _make_user("admin1234", "admin", "Clinical Admin"),
    "analyst": _make_user("analyst1234", "analyst", "Clinical Analyst"),
    "viewer": _make_user("viewer1234", "viewer", "Read-Only Viewer"),
}

def _load_json_file(path: Path, default):
    try:
        if path.exists():
            with path.open("r", encoding="utf-8") as f:
                return json.load(f)
    except Exception:
        pass
    return default

def _write_json_file(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    tmp.replace(path)

def _normalize_user_record(username: str, user: dict) -> dict:
    now = datetime.utcnow().isoformat()
    user = dict(user or {})
    # Migration support for old sha256-only stores, if ever present.
    if "hash" in user and ("salt" not in user or "dk" not in user):
        fallback_password = f"{username}1234" if username in {"admin", "analyst", "viewer"} else secrets.token_urlsafe(12)
        salt, dk = _hash_password(fallback_password)
        user["salt"], user["dk"] = salt, dk
        user.pop("hash", None)
    user.setdefault("role", "viewer")
    user.setdefault("display", username)
    user.setdefault("is_active", True)
    user.setdefault("created_at", now)
    user.setdefault("updated_at", now)
    user.setdefault("last_login", None)
    return user

def _load_users() -> Dict[str, dict]:
    loaded = _load_json_file(USER_STORE_PATH, None)
    raw = loaded if isinstance(loaded, dict) and loaded else DEFAULT_USERS
    users = {str(k).strip().lower(): _normalize_user_record(str(k).strip().lower(), v) for k, v in raw.items()}
    if "admin" not in users:
        users["admin"] = DEFAULT_USERS["admin"]
    return users

def _save_users() -> None:
    _write_json_file(USER_STORE_PATH, USERS)

def _audit(actor: str, action: str, target: str = "", details: Optional[dict] = None) -> None:
    log = _load_json_file(AUDIT_LOG_PATH, [])
    if not isinstance(log, list):
        log = []
    log.append({
        "timestamp": datetime.utcnow().isoformat(),
        "actor": actor,
        "action": action,
        "target": target,
        "details": details or {},
    })
    _write_json_file(AUDIT_LOG_PATH, log[-250:])

def _public_user(username: str, user: dict) -> dict:
    return {
        "username": username,
        "role": user.get("role", "viewer"),
        "display": user.get("display", username),
        "is_active": bool(user.get("is_active", True)),
        "created_at": user.get("created_at"),
        "updated_at": user.get("updated_at"),
        "last_login": user.get("last_login"),
    }

def _active_admin_count(exclude_username: Optional[str] = None) -> int:
    return sum(
        1 for uname, user in USERS.items()
        if uname != exclude_username and user.get("role") == "admin" and user.get("is_active", True)
    )

def _invalidate_user_sessions(username: str, keep_current_token: Optional[str] = None) -> None:
    for token in list(SESSIONS.keys()):
        if SESSIONS[token].get("username") == username and token != keep_current_token:
            del SESSIONS[token]

USERS: Dict[str, dict] = _load_users()
SESSIONS: Dict[str, dict] = {}

def _purge_expired() -> None:
    now = datetime.utcnow()
    for token in [t for t, sess in SESSIONS.items() if sess.get("expires") and sess["expires"] <= now]:
        del SESSIONS[token]

def _get_session(session_id: Optional[str] = Depends(cookie_scheme)):
    _purge_expired()
    if not session_id or session_id not in SESSIONS:
        raise HTTPException(status_code=401, detail="Not authenticated")
    sess = SESSIONS[session_id]
    if sess["expires"] <= datetime.utcnow():
        del SESSIONS[session_id]
        raise HTTPException(status_code=401, detail="Session expired")
    return sess

def _require_role(minimum_role: str):
    def _check(session: dict = Depends(_get_session)) -> dict:
        role = session.get("role", "viewer")
        if ROLE_RANK.get(role, 0) < ROLE_RANK[minimum_role]:
            raise HTTPException(status_code=403, detail=f"Access denied: minimum role is {minimum_role}")
        return session
    return _check

AnyAuth = Depends(_get_session)
AnalystUp = Depends(_require_role("analyst"))
AdminOnly = Depends(_require_role("admin"))

@app.post("/api/login")
def login(data: dict, request: Request):
    username = (data.get("username") or "").strip().lower()
    password = data.get("password") or ""
    user = USERS.get(username)
    if not user or not _verify_password(password, user["salt"], user["dk"]):
        raise HTTPException(status_code=401, detail="Invalid credentials")
    if not user.get("is_active", True):
        raise HTTPException(status_code=403, detail="Account is disabled")

    token = secrets.token_urlsafe(32)
    user["last_login"] = datetime.utcnow().isoformat()
    user["updated_at"] = user.get("updated_at") or user["last_login"]
    _save_users()
    _audit(username, "login", username)
    SESSIONS[token] = {
        "username": username,
        "role": user["role"],
        "display": user["display"],
        "expires": datetime.utcnow() + timedelta(hours=SESSION_TTL_HOURS),
    }
    resp = JSONResponse({"username": username, "role": user["role"], "display": user["display"]})
    forwarded_proto = request.headers.get("x-forwarded-proto", "")
    cookie_secure = request.url.scheme == "https" or forwarded_proto == "https"
    resp.set_cookie(key="session_id", value=token, httponly=True, secure=cookie_secure, samesite="strict", max_age=SESSION_TTL_HOURS * 3600)
    return resp

@app.post("/api/logout")
def logout(session: dict = Depends(_get_session), session_id: Optional[str] = Depends(cookie_scheme)):
    if session_id and session_id in SESSIONS:
        del SESSIONS[session_id]
    resp = JSONResponse({"ok": True})
    resp.delete_cookie("session_id")
    return resp

@app.get("/api/me")
def me(session: dict = AnyAuth):
    return {"username": session["username"], "role": session["role"], "display": session["display"]}

@app.get("/api/admin/users")
def admin_list_users(session: dict = AdminOnly):
    users = [_public_user(username, user) for username, user in USERS.items()]
    users.sort(key=lambda x: (x["role"] != "admin", x["username"]))
    return {"users": users}

@app.post("/api/admin/users")
def admin_create_user(data: dict, session: dict = AdminOnly):
    username = (data.get("username") or "").strip().lower()
    display = (data.get("display") or username).strip()
    role = (data.get("role") or "viewer").strip().lower()
    password = str(data.get("password") or "")
    if not re.fullmatch(r"[a-z0-9_.-]{3,32}", username):
        raise HTTPException(status_code=400, detail="Username must be 3-32 characters: a-z, 0-9, dot, underscore, or hyphen")
    if username in USERS:
        raise HTTPException(status_code=409, detail="Username already exists")
    if role not in ROLE_RANK:
        raise HTTPException(status_code=400, detail="Invalid role")
    if len(password) < 8:
        raise HTTPException(status_code=400, detail="Password must be at least 8 characters")
    USERS[username] = _make_user(password, role, display)
    USERS[username]["is_active"] = bool(data.get("is_active", True))
    _save_users()
    _audit(session["username"], "create_user", username, {"role": role, "is_active": USERS[username]["is_active"]})
    return {"ok": True, "user": _public_user(username, USERS[username])}

@app.post("/api/admin/users/{username}")
def admin_update_user(username: str, data: dict, session: dict = AdminOnly, session_id: Optional[str] = Depends(cookie_scheme)):
    username = username.strip().lower()
    if username not in USERS:
        raise HTTPException(status_code=404, detail="User not found")
    current = USERS[username]
    role = (data.get("role") or current.get("role") or "viewer").strip().lower()
    if role not in ROLE_RANK:
        raise HTTPException(status_code=400, detail="Invalid role")
    display = (data.get("display") or current.get("display") or username).strip()
    is_active = bool(data.get("is_active", current.get("is_active", True)))
    if username == session["username"] and (role != "admin" or not is_active):
        raise HTTPException(status_code=400, detail="You cannot remove your own admin access or disable your own account")
    if current.get("role") == "admin" and current.get("is_active", True) and (role != "admin" or not is_active):
        if _active_admin_count(exclude_username=username) == 0:
            raise HTTPException(status_code=400, detail="At least one active admin is required")
    current["role"] = role
    current["display"] = display
    current["is_active"] = is_active
    current["updated_at"] = datetime.utcnow().isoformat()
    password = data.get("password")
    if password:
        password = str(password)
        if len(password) < 8:
            raise HTTPException(status_code=400, detail="Password must be at least 8 characters")
        current["salt"], current["dk"] = _hash_password(password)
        _invalidate_user_sessions(username, keep_current_token=session_id if username == session["username"] else None)
    if not is_active:
        _invalidate_user_sessions(username, keep_current_token=session_id if username == session["username"] else None)
    for sess in SESSIONS.values():
        if sess.get("username") == username:
            sess["role"] = role
            sess["display"] = display
    _save_users()
    _audit(session["username"], "update_user", username, {"role": role, "is_active": is_active, "password_changed": bool(password)})
    return {"ok": True, "user": _public_user(username, current)}

@app.delete("/api/admin/users/{username}")
def admin_delete_user(username: str, session: dict = AdminOnly):
    username = username.strip().lower()
    if username not in USERS:
        raise HTTPException(status_code=404, detail="User not found")
    if username == session["username"]:
        raise HTTPException(status_code=400, detail="You cannot delete your own account")
    user = USERS[username]
    if user.get("role") == "admin" and user.get("is_active", True) and _active_admin_count(exclude_username=username) == 0:
        raise HTTPException(status_code=400, detail="At least one active admin is required")
    del USERS[username]
    _invalidate_user_sessions(username)
    _save_users()
    _audit(session["username"], "delete_user", username)
    return {"ok": True}

@app.get("/api/admin/audit")
def admin_audit_log(session: dict = AdminOnly):
    log = _load_json_file(AUDIT_LOG_PATH, [])
    if not isinstance(log, list):
        log = []
    return {"events": list(reversed(log[-50:]))}

# --- Data Engine ---

def _fs():
    if USE_LOCAL_DATA: return None
    import gcsfs
    return gcsfs.GCSFileSystem()

def _load_csv(path: str) -> pd.DataFrame:
    try:
        if USE_LOCAL_DATA:
            p = DATA_DIR / path.replace("/", os.sep)
            if not p.exists():
                return pd.DataFrame()
            df = pd.read_csv(p)
        else:
            full_path = f"gs://{BUCKET}/{path}"
            df = pd.read_csv(full_path)
        
        df.columns = [c.strip() for c in df.columns]
        if "pid" not in df.columns:
            if "sample" in df.columns: df["pid"] = df["sample"]
            elif "PID" in df.columns: df["pid"] = df["PID"]
            
        for c in ["patient_key", "sample_id", "pid", "sample"]:
            if c in df.columns:
                df[c] = df[c].astype(str).str.strip()
        
        if "timepoint" in df.columns:
            df["timepoint"] = pd.to_numeric(df["timepoint"], errors="coerce").astype("Int64")
        return df
    except Exception:
        return pd.DataFrame()

def _master():
    return _load_csv(MASTER_CSV)

def _first_nonnull(s):
    v = s.dropna()
    return v.iloc[0] if not v.empty else None

def _clean_json(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _clean_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_clean_json(v) for v in obj]
    elif isinstance(obj, (float, np.float64, np.float32)):
        if np.isnan(obj) or np.isinf(obj):
            return None
        return float(obj)
    return obj

def _json_records(df: pd.DataFrame):
    if df.empty: return []
    return _clean_json(df.to_dict(orient="records"))

def _list_reports():
    cache_key = "reports_cache"
    if cache_key in _cache: return _cache[cache_key]
    if USE_LOCAL_DATA:
        lp = DATA_DIR / "reports"
        if lp.exists():
             names = [f.name for f in lp.glob("*.pdf")]
             _cache[cache_key] = names
             return names
        return []
    fs = _fs()
    if fs is None: return []
    try:
        files = fs.ls(f"{BUCKET}/{REPORTS_DIR}")
        names = [Path(f).name for f in files]
        _cache[cache_key] = names
        return names
    except Exception: return []

@app.get("/api/overview")
def overview(session: dict = Depends(_get_session)):
    m = _master()
    if m.empty or "pid" not in m.columns:
        return JSONResponse({"patients": [], "stats": {}})
    
    m = m[m["pid"].astype(str).str.startswith("2LB")].copy()
    
    base_info = m.dropna(subset=["pid"]).groupby("pid")["timepoint"].apply(lambda s: sorted([int(x) for x in s.dropna().unique()])).reset_index(name="timepoints")
    base_info["n_timepoints"] = base_info["timepoints"].apply(len)
    base_info = base_info.rename(columns={"pid": "sample"}) 

    clin_fields = ["age", "Gender", "Treatment", "TreatmentGroup", "CT", "pfs_days", "egfr_status"]
    for f in clin_fields:
        if f in m.columns:
            agg = m.groupby("pid")[f].apply(_first_nonnull).reset_index(name=f)
            base_info = base_info.merge(agg, left_on="sample", right_on="pid", how="left").drop(columns=["pid"])
            
    if "tumor_fraction" in m.columns:
        tf_agg = m.groupby("pid")["tumor_fraction"].mean().reset_index(name="mean_tf")
        base_info = base_info.merge(tf_agg, left_on="sample", right_on="pid", how="left").drop(columns=["pid"])
    
    try:
        ml = _load_csv(ML_REGISTRY_CSV)
        if not ml.empty:
            ml = ml[ml["pid"].astype(str).str.startswith("2LB")].copy()
            ml_latest = ml.sort_values(["pid", "timepoint"]).groupby("pid").tail(1)
            ml_latest["pred"] = ml_latest["response_prob"].apply(lambda p: "PR" if (p or 0) > 0.5 else "non-PR")
            ml_latest["actual"] = ml_latest["CT"].apply(lambda ct: "PR" if str(ct).strip().upper() == "PR" else "non-PR")
            
            cols = ["cp_score", "response_prob", "pred", "actual", "pfs_days", "age", "Gender", "egfr_status"]
            valid_ml = [c for c in cols if c in ml_latest.columns]
            ml_map = ml_latest.set_index("pid")[valid_ml]
            base_info = base_info.merge(ml_map, left_on="sample", right_index=True, how="left", suffixes=("", "_ml"))
            
            ml_counts = ml.groupby("pid").size().reset_index(name="n_ml_timepoints")
            base_info = base_info.merge(ml_counts, left_on="sample", right_on="pid", how="left").drop(columns=["pid"])
            base_info["n_ml_timepoints"] = base_info["n_ml_timepoints"].fillna(0).astype(int)
    except Exception: base_info["n_ml_timepoints"] = 0

    if "TreatmentGroup" in base_info.columns:
        def map_tx(x):
            xs = str(x).strip()
            if xs.lower() in ["control/unknown", "control", "none"]: return "null"
            if xs in ["TKI", "Targeted (Other)"]: return "Targeted Therapy"
            if xs == "Chemo-IO": return "Chemotherapy"
            return xs
        base_info["TreatmentGroup"] = base_info["TreatmentGroup"].apply(map_tx)
    
    for f in ["age", "pfs_days", "Gender", "egfr_status"]:
        ml_f = f + "_ml"
        if ml_f in base_info.columns:
            base_info[f] = base_info[f].fillna(base_info[ml_f])
            base_info.drop(columns=[ml_f], inplace=True)
    
    base_info["pfs_val"] = base_info.get("pfs_days")
    
    reports = _list_reports()
    if "patient_key" in m.columns:
        m["_has_pdf"] = m["patient_key"].apply(lambda pk: f"{pk}_genomeWide.pdf" in reports)
        pdf_agg = m.groupby("pid")["_has_pdf"].sum().reset_index(name="has_pdf_count")
        base_info = base_info.merge(pdf_agg, left_on="sample", right_on="pid", how="left").drop(columns=["pid"])

    stats = {"bucket": BUCKET}
    if not base_info.empty:
        if "CT" in base_info.columns: stats["ct_dist"] = base_info["CT"].value_counts().to_dict()
        if "TreatmentGroup" in base_info.columns: stats["tx_dist"] = base_info["TreatmentGroup"].value_counts().to_dict()
        if "Gender" in base_info.columns:
            stats["male_count"] = int((base_info["Gender"] == "Male").sum())
            stats["female_count"] = int((base_info["Gender"] == "Female").sum())
        if "mean_tf" in base_info.columns:
            tfs = pd.to_numeric(base_info["mean_tf"], errors="coerce").dropna()
            if not tfs.empty:
                stats["median_tf"] = float(tfs.median())
                stats["tf_dist"] = {str(k): int(v) for k, v in pd.cut(tfs, bins=[0, 0.01, 0.05, 0.1, 0.5, 1.0]).value_counts().sort_index().items()}
        if "response_prob" in base_info.columns:
            rps = pd.to_numeric(base_info["response_prob"], errors="coerce").dropna()
            if not rps.empty: stats["mean_rp"] = float(rps.mean())
        if all(x in base_info.columns for x in ["CT", "cp_score"]):
            stats["urgent_count"] = int(((base_info["CT"] == "PD") | (base_info["cp_score"] > 0.8)).sum())
            stats["discordant_count"] = int(((base_info["CT"] == "PR") & (base_info["cp_score"] > 0.8)).sum())
        
    return JSONResponse({
        "patients": _json_records(base_info.sort_values("sample")),
        "stats": _clean_json(stats),
        "total_samples": int(m.shape[0]),
        "total_patients": int(base_info.shape[0])
    })

@app.get("/api/patient/{sample_id}")
def patient_timepoints(sample_id: str, session: dict = Depends(_get_session)):
    if not sample_id.startswith("2LB"): return JSONResponse([])
    m = _master()
    sub = m[m["pid"].astype(str) == sample_id].copy()
    if sub.empty: return JSONResponse([])
    reports = _list_reports()
    if "patient_key" in sub.columns: sub["has_pdf"] = sub["patient_key"].apply(lambda pk: f"{pk}_genomeWide.pdf" in reports)
    try:
        ml = _load_csv(ML_REGISTRY_CSV)
        if not ml.empty:
            ml["pred"] = ml["response_prob"].apply(lambda p: "PR" if (p or 0) > 0.5 else "non-PR")
            ml["actual"] = ml["CT"].apply(lambda ct: "PR" if str(ct).strip().upper() == "PR" else "non-PR")
            p_cols = ["sample_id", "cp_score", "response_prob", "pred", "actual", "pfs_days", "age", "egfr_status"]
            v_p_cols = [c for c in p_cols if c in ml.columns]
            sub = sub.merge(ml[v_p_cols], on="sample_id", how="left", suffixes=("", "_ml"))
    except Exception: pass
    for f in ["age", "pfs_days", "egfr_status"]:
        if f in sub.columns: sub[f] = sub[f].fillna(method="ffill").fillna(method="bfill")
    sub["pfs_val"] = sub.get("pfs_days_ml", sub.get("pfs_days"))
    return JSONResponse(_json_records(sub))

@app.get("/api/features/{patient_key}/fraglen")
def fraglen(patient_key: str, session: dict = Depends(_get_session)):
    df = _load_csv(FRAG_BINS_CSV)
    if df.empty or "patient_key" not in df.columns: return JSONResponse([])
    
    row = df[df["patient_key"].astype(str) == patient_key]
    is_fallback = False
    if row.empty:
        # Fallback: Extract PID prefix (e.g., 2LB-001) and find latest available for that patient
        pid_prefix = patient_key.split('P-')[0] if 'P-' in patient_key else patient_key
        patient_rows = df[df["patient_key"].astype(str).str.startswith(pid_prefix)]
        if not patient_rows.empty:
            row = patient_rows.sort_values("timepoint", ascending=False).head(1)
            is_fallback = True
            
    if row.empty: return JSONResponse([])
    
    row_data = row.iloc[0]; result = []
    for c in row_data.index:
        c_str = str(c)
        if c_str.startswith("bin_"):
            parts = c_str.split("_")
            if len(parts) >= 3:
                try:
                    x_val = (int(parts[1]) + int(parts[2])) / 2
                    result.append({"x": x_val, "y": float(row_data[c])})
                except: pass
        elif c_str.isdigit(): 
            try: result.append({"x": int(c_str), "y": float(row_data[c])})
            except: pass
            
    sorted_res = sorted(result, key=lambda r: r["x"])
    return JSONResponse(_clean_json({"data": sorted_res, "fallback": is_fallback, "source_key": str(row_data.get("patient_key", ""))}))

@app.get("/api/features/{patient_key}/endmotif")
def endmotif(patient_key: str, topk: int = 12, session: dict = Depends(_get_session)):
    em = _load_csv(ENDMOTIF_CSV)
    if em.empty or "patient_key" not in em.columns: return JSONResponse([])
    
    row = em[em["patient_key"].astype(str) == patient_key]
    is_fallback = False
    if row.empty:
        pid_prefix = patient_key.split('P-')[0] if 'P-' in patient_key else patient_key
        patient_rows = em[em["patient_key"].astype(str).str.startswith(pid_prefix)]
        if not patient_rows.empty:
            row = patient_rows.sort_values("timepoint", ascending=False).head(1)
            is_fallback = True

    if row.empty: return JSONResponse([])
    
    row_data = row.iloc[0]; vals = []
    for c in row_data.index:
        if re.fullmatch(r"[ACGT]{4}", str(c)) and not pd.isna(row_data[c]): 
            vals.append({"motif": str(c), "value": float(row_data[c])})
    
    vals.sort(key=lambda x: x["value"], reverse=True)
    return JSONResponse(_clean_json({"data": vals[:topk], "fallback": is_fallback, "source_key": str(row_data.get("patient_key", ""))}))

@app.get("/api/report/cnv/{sample_id}")
@app.get("/api/patient/{patient_key}/pdf")
def get_pdf(patient_key: Optional[str] = None, sample_id: Optional[str] = None, session: dict = AnalystUp):
    pk = sample_id or patient_key
    print(f"DEBUG: get_pdf requested for: {pk}")
    if not pk: raise HTTPException(status_code=400, detail="Missing ID")
    fname = f"{pk}_genomeWide.pdf"
    if USE_LOCAL_DATA:
        lp = DATA_DIR / "reports" / fname
        if lp.exists(): return FileResponse(lp, media_type="application/pdf", filename=fname)
        raise HTTPException(status_code=404, detail="CNA PDF not found locally")
    fs = _fs()
    if fs:
        try:
            full_p = f"{BUCKET}/{REPORTS_DIR}/{fname}"
            with fs.open(full_p, "rb") as f: content = f.read()
            return Response(content=content, media_type="application/pdf", headers={"Content-Disposition": f'inline; filename="{fname}"'})
        except Exception: pass
    raise HTTPException(status_code=404, detail="CNA PDF not found")

EXPORT_VERSION = "v4.5"
DISCLAIMER = "DISCLAIMER: ML predictions are restricted to the monitoring cohort. Research use only."

def _add_patient_to_zip(z: zipfile.ZipFile, pid: str, root_folder: str = ""):
    print(f"DEBUG: Bundling patient {pid} into {root_folder}")
    m = _master(); pv = m[m["pid"] == pid].sort_values("timepoint")
    if pv.empty: return
    ml = _load_csv(ML_REGISTRY_CSV); pml = ml[ml["pid"] == pid].sort_values("timepoint")
    reports = _list_reports(); summary_df = pv.merge(pml, on="sample_id", how="left", suffixes=("", "_ml"))
    summary_df["has_pdf"] = summary_df["sample_id"].apply(lambda sid: f"{sid}_genomeWide.pdf" in reports)
    
    # Fragmentomics & Motifs for the ZIP
    frag_all = _load_csv(FRAG_BINS_CSV); p_frag = frag_all[frag_all["patient_key"].isin(pv["patient_key"])]
    motif_all = _load_csv(ENDMOTIF_CSV); p_motif = motif_all[motif_all["patient_key"].isin(pv["patient_key"])]

    profile = {
        "pid": pid, 
        "ct_label": _first_nonnull(pv.get("CT", pd.Series([None]))), 
        "age": _first_nonnull(summary_df.get("age", pd.Series([None]))), 
        "version": EXPORT_VERSION,
        "export_timestamp": datetime.utcnow().isoformat()
    }
    
    prefix = f"{root_folder}/" if root_folder else ""
    z.writestr(f"{prefix}longitudinal_summary.csv", summary_df.to_csv(index=False))
    z.writestr(f"{prefix}patient_profile.json", pd.Series(profile).to_json())
    z.writestr(f"{prefix}ml_evidence.csv", pml.to_csv(index=False))
    z.writestr(f"{prefix}fragmentomics.csv", p_frag.to_csv(index=False))
    z.writestr(f"{prefix}end_motifs.csv", p_motif.to_csv(index=False))
    z.writestr(f"{prefix}README.txt", f"Patient ID: {pid}\nExported: {profile['export_timestamp']}\n\n{DISCLAIMER}")
    
    for _, row in summary_df[summary_df["has_pdf"]].iterrows():
        sid = row["sample_id"]; fn = f"{sid}_genomeWide.pdf"
        try:
            if USE_LOCAL_DATA:
                lp = DATA_DIR / "reports" / fn
                if lp.exists(): z.writestr(f"{prefix}CNA_Reports/{fn}", lp.read_bytes())
            else:
                fs = _fs()
                if fs:
                    with fs.open(f"{BUCKET}/{REPORTS_DIR}/{fn}", "rb") as f: z.writestr(f"{prefix}CNA_Reports/{fn}", f.read())
        except: pass

@app.get("/api/export/patient/{pid}")
def export_patient(pid: str, session: dict = AnalystUp):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as z:
        _add_patient_to_zip(z, pid, root_folder="") # Flattened for single bundle
    return Response(content=buf.getvalue(), media_type="application/zip", headers={"Content-Disposition": f'attachment; filename="{pid}_bundle.zip"'})

@app.post("/api/export/selected")
def export_selected(data: dict, session: dict = AdminOnly):
    pids = data.get("pids", [])
    if not pids: raise HTTPException(400, "No patients selected")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for pid in pids: _add_patient_to_zip(z, pid, root_folder=pid)
        z.writestr("selected_manifest.json", pd.Series({"pids": pids, "timestamp": datetime.utcnow().isoformat()}).to_json())
    return Response(content=buf.getvalue(), media_type="application/zip", headers={"Content-Disposition": 'attachment; filename="selected_patient_bundles.zip"'})

@app.get("/api/export/all")
@app.post("/api/export/all")
def export_all(session: dict = AdminOnly):
    m = _master()
    pids = m[m["pid"].astype(str).str.startswith("2LB")]["pid"].unique().tolist()
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for pid in pids: _add_patient_to_zip(z, pid, root_folder=pid)
        z.writestr("cohort_manifest.json", pd.Series({"cohort": "2LB", "count": len(pids), "timestamp": datetime.utcnow().isoformat()}).to_json())
    return Response(content=buf.getvalue(), media_type="application/zip", headers={"Content-Disposition": 'attachment; filename="all_cohort_patient_bundles.zip"'})

FRONTEND = os.path.join(os.path.dirname(__file__), "..", "frontend")
@app.get("/{full_path:path}")
def spa(full_path: str):
    static = os.path.join(FRONTEND, full_path)
    if os.path.isfile(static): return FileResponse(static)
    return FileResponse(os.path.join(FRONTEND, "index.html"))
