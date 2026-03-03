# nenufar_sb_scan.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any
import json
import glob
import json
import glob
import re

def _get_table():
    try:
        from casacore.tables import table
        return table, "casacore"
    except Exception:
        try:
            from pyrap.tables import table
            return table, "pyrap"
        except Exception as e:
            raise RuntimeError("Neither python-casacore nor pyrap is available.") from e

def _ms_center_freq_mhz(ms_path: Path, table) -> float:
    with table(str(ms_path) + "::SPECTRAL_WINDOW", ack=False) as t:
        f = t.getcol("CHAN_FREQ")  # Hz
        return float(f.mean() / 1e6)

@dataclass
class ScanResult:
    base_dir: str
    event_dir: str
    l1_dir: str
    freq_min_mhz: Optional[float]
    freq_max_mhz: Optional[float]
    backend: str
    n_total: int
    n_ok: int
    n_bad: int
    csv_all: str
    csv_sel: str
    json_sel: str
    selected_sbs: List[str]

def scan_sb_freq(
    l1_dir: str | Path,
    workdir: str | Path,
    freq_range_mhz: Optional[Tuple[float, float]] = (70.0, 80.0),
    pattern: str = "SB*.MS",
) -> Tuple["Any", ScanResult]:
    """
    Scan all SB*.MS under l1_dir, compute center frequency, save CSV/JSON.

    freq_range_mhz:
      - None -> select all
      - (fmin, fmax) -> filter inclusive
    """
    import pandas as pd

    l1_dir = Path(l1_dir)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    table, backend = _get_table()

    sb_paths = sorted(l1_dir.glob(pattern))
    rows: List[Dict[str, Any]] = []
    for ms in sb_paths:
        try:
            ctr = _ms_center_freq_mhz(ms, table)
            rows.append({"sb": ms.name, "ms": str(ms), "ctr_mhz": ctr})
        except Exception as e:
            rows.append({"sb": ms.name, "ms": str(ms), "ctr_mhz": None, "error": str(e)})

    rows_ok = [r for r in rows if r.get("ctr_mhz") is not None]
    rows_ok = sorted(rows_ok, key=lambda r: r["ctr_mhz"])
    rows_bad = [r for r in rows if r.get("ctr_mhz") is None]

    df = pd.DataFrame(rows_ok)

    csv_all = workdir / "sb_freq_table_all.csv"
    df.to_csv(csv_all, index=False)

    if freq_range_mhz is None:
        df_sel = df.copy()
    else:
        fmin, fmax = freq_range_mhz
        df_sel = df[(df["ctr_mhz"] >= fmin) & (df["ctr_mhz"] <= fmax)].copy()

    if freq_range_mhz is None:
        csv_sel = workdir / "sb_freq_table_selected_ALL.csv"
    else:
        csv_sel = workdir / f"sb_freq_table_selected_{freq_range_mhz[0]:.0f}-{freq_range_mhz[1]:.0f}MHz.csv"
    df_sel.to_csv(csv_sel, index=False)

    selected_ms = df_sel["ms"].tolist()
    selected_sb = df_sel["sb"].tolist()

    json_sel = workdir / "selected_sb_list.json"
    json.dump(
        {
            "freq_range_mhz": None if freq_range_mhz is None else [float(freq_range_mhz[0]), float(freq_range_mhz[1])],
            "selected_ms": selected_ms,
            "selected_sb": selected_sb,
        },
        open(json_sel, "w"),
        indent=2,
    )

    meta = ScanResult(
        base_dir=str(l1_dir.parent.parent.parent.parent.parent) if len(l1_dir.parts) > 6 else "",
        event_dir=str(l1_dir.parent),
        l1_dir=str(l1_dir),
        freq_min_mhz=None if freq_range_mhz is None else float(freq_range_mhz[0]),
        freq_max_mhz=None if freq_range_mhz is None else float(freq_range_mhz[1]),
        backend=backend,
        n_total=len(sb_paths),
        n_ok=len(rows_ok),
        n_bad=len(rows_bad),
        csv_all=str(csv_all),
        csv_sel=str(csv_sel),
        json_sel=str(json_sel),
        selected_sbs=selected_sb,
    )
    return df_sel, meta

def list_sun_tracking_events(base_dir: str | Path) -> List[str]:
    """
    Find SUN_TRACKING directories under:
      base_dir/YYYY/MM/
    Returns sorted full paths to event dirs.
    """
    base_dir = Path(base_dir)
    events = sorted(glob.glob(str(base_dir / "*" / "*" / "*SUN_TRACKING")))
    return events

def infer_l1_dir(event_dir: str | Path) -> Path:
    event_dir = Path(event_dir)
    l1 = event_dir / "L1"
    if not l1.exists():
        raise FileNotFoundError(f"L1 not found: {l1}")
    return l1

def scan_by_ymd(base, workroot, year, month, day, freq_range_mhz=None):
    """
    Find SUN_TRACKING event on a given date, then scan its L1.
    If multiple SUN_TRACKING events exist that day, pick the longest-duration one.
    """
    from pathlib import Path
    import glob
    import re

    base = Path(base)
    workroot = Path(workroot)

    yy = f"{int(year):04d}"
    mm = f"{int(month):02d}"
    dd = f"{int(day):02d}"

    candidates = sorted(glob.glob(str(base / yy / mm / f"{yy}{mm}{dd}_*SUN_TRACKING")))
    if not candidates:
        raise FileNotFoundError(f"No SUN_TRACKING found for {yy}-{mm}-{dd} under {base/yy/mm}")

    # choose longest (roughly) using HHMMSS difference in folder name
    def _dur_key(p):
        name = Path(p).name
        match = re.search(
            rf"{yy}{mm}{dd}_(\d{{6}})_{yy}{mm}{dd}_(\d{{6}})_SUN_TRACKING",
            name
        )
        if not match:
            return 0
        t1, t2 = match.group(1), match.group(2)
        return int(t2) - int(t1)

    best = sorted(candidates, key=_dur_key, reverse=True)[0]
    l1 = infer_l1_dir(best)
    workdir = Path(workroot) / Path(best).name
    return scan_sb_freq(l1, workdir=workdir, freq_range_mhz=freq_range_mhz)

def pick_event_dir(base, year, month, day, suffix):
    """
    suffix: "SUN_TRACKING" or "CAS_A_TRACKING"
    Pick the longest-duration event folder for that day.
    """
    from pathlib import Path
    import glob, re

    base = Path(base)
    yy = f"{int(year):04d}"
    mm = f"{int(month):02d}"
    dd = f"{int(day):02d}"

    candidates = sorted(glob.glob(str(base / yy / mm / f"{yy}{mm}{dd}_*_{suffix}")))
    if not candidates:
        raise FileNotFoundError(f"No {suffix} found for {yy}-{mm}-{dd} under {base/yy/mm}")

    def _dur_key(p):
        name = Path(p).name
        match = re.search(rf"{yy}{mm}{dd}_(\d{{6}})_{yy}{mm}{dd}_(\d{{6}})_{suffix}$", name)
        if not match:
            return 0
        return int(match.group(2)) - int(match.group(1))

    return sorted(candidates, key=_dur_key, reverse=True)[0]


def scan_sun_and_casa_by_ymd(base, workroot, year, month, day, freq_range_mhz=None, casa_event_dir=None):
    """
    Output:
      df with: sb, ctr_mhz, sun_ms,
               casa_pre_ms, casa_post_ms,
               casa_chosen_ms, casa_chosen_tag
    Also meta includes all CASA candidates + chosen one.
    """
    from pathlib import Path
    import json

    # 1) pick SUN event
    sun_event = pick_event_dir(base, year, month, day, "SUN_TRACKING")

    # 2) list ALL CASA candidates that day
    casa_candidates = list_event_dirs_by_date(base, year, month, day, "CAS_A_TRACKING")
    if not casa_candidates:
        raise FileNotFoundError("No CAS_A_TRACKING candidates found for this date.")

    # 3) split into pre/post relative to SUN start
    casa_pre, casa_post = split_casa_candidates_relative_to_sun(sun_event, casa_candidates)

    # 4) choose CASA event
    if casa_event_dir is None:
        # default: closest before; if none, closest after
        casa_event = pick_closest_calibrator(sun_event, casa_candidates)
    else:
        casa_event = casa_event_dir

    # Determine chosen tag (pre/post/unknown)
    casa_chosen_tag = "unknown"
    if str(casa_event) in [str(x) for x in casa_pre]:
        casa_chosen_tag = "pre"
    elif str(casa_event) in [str(x) for x in casa_post]:
        casa_chosen_tag = "post"

    # 5) resolve L1 dirs (chosen + optional pre/post)
    sun_l1 = infer_l1_dir(sun_event)
    casa_l1_chosen = infer_l1_dir(casa_event)

    casa_l1_pre = infer_l1_dir(casa_pre[-1]) if casa_pre else None   # usually only one
    casa_l1_post = infer_l1_dir(casa_post[0]) if casa_post else None # usually only one

    # 6) scan SUN SBs
    df_sel, _meta = scan_sb_freq(sun_l1, workdir=Path(workroot)/Path(sun_event).name, freq_range_mhz=freq_range_mhz)
    df_sel = df_sel.copy()
    df_sel["sun_ms"] = df_sel["ms"]

    # 7) build CASA paths
    def _mk(msdir, sb):
        return str(Path(msdir)/sb) if msdir is not None else None

    df_sel["casa_pre_ms"]  = df_sel["sb"].apply(lambda sb: _mk(casa_l1_pre, sb))
    df_sel["casa_post_ms"] = df_sel["sb"].apply(lambda sb: _mk(casa_l1_post, sb))
    df_sel["casa_chosen_ms"] = df_sel["sb"].apply(lambda sb: str(Path(casa_l1_chosen)/sb))

    # existence checks (optional but useful)
    df_sel["casa_chosen_exists"] = df_sel["casa_chosen_ms"].apply(lambda p: Path(p).exists())
    df_sel = df_sel[df_sel["casa_chosen_exists"]].reset_index(drop=True)

    # 8) write JSON for downstream DP3
    outdir = Path(workroot) / Path(sun_event).name
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / "selected_sb_pair_list.json"

    payload = {
        "date": f"{int(year):04d}-{int(month):02d}-{int(day):02d}",
        "freq_range_mhz": None if freq_range_mhz is None else [float(freq_range_mhz[0]), float(freq_range_mhz[1])],
        "sun_event": str(sun_event),
        "casa_candidates": [str(x) for x in casa_candidates],
        "casa_pre": [str(x) for x in casa_pre],
        "casa_post": [str(x) for x in casa_post],
        "casa_event_chosen": str(casa_event),
        "casa_chosen_tag": casa_chosen_tag,
        "sun_l1": str(sun_l1),
        "casa_l1_chosen": str(casa_l1_chosen),
        "selected_sb": df_sel["sb"].tolist(),
        "ctr_mhz": df_sel["ctr_mhz"].tolist(),
        "sun_ms": df_sel["sun_ms"].tolist(),
        "casa_ms": df_sel["casa_chosen_ms"].tolist(),  # downstream DP3 uses this
    }
    json.dump(payload, open(json_path, "w"), indent=2)

    meta_out = {
        "sun_event": str(sun_event),
        "casa_candidates": [str(x) for x in casa_candidates],
        "casa_event_chosen": str(casa_event),
        "casa_chosen_tag": casa_chosen_tag,
        "n_sb_selected": int(len(df_sel)),
        "pair_json": str(json_path),
        "workdir": str(outdir),
    }

    return df_sel[["sb","ctr_mhz","sun_ms","casa_pre_ms","casa_post_ms","casa_chosen_ms"]], meta_out

def list_event_dirs_by_date(base, year, month, day, suffix):
    """
    Return all event dirs matching yyyymmdd_*_SUFFIX under base/YYYY/MM.
    suffix examples: "SUN_TRACKING", "CAS_A_TRACKING"
    """
    base = Path(base)
    yy = f"{int(year):04d}"
    mm = f"{int(month):02d}"
    dd = f"{int(day):02d}"
    return sorted(glob.glob(str(base / yy / mm / f"{yy}{mm}{dd}_*_{suffix}")))

def pick_closest_calibrator(sun_event_dir, casa_event_dirs):
    """
    Default strategy: choose calibrator whose end time is closest BEFORE Sun start;
    if none before, choose the closest AFTER.
    """
    sun_name = Path(sun_event_dir).name
    m = re.search(r"_(\d{6})_\d{8}_\d{6}_SUN_TRACKING$", sun_name)
    # fallback: parse by splitting
    sun_start = None
    try:
        # 20240310_101000_20240310_135000_SUN_TRACKING
        sun_start = int(sun_name.split("_")[1])
    except Exception:
        sun_start = None

    def parse_times(ev):
        name = Path(ev).name
        try:
            t1 = int(name.split("_")[1])
            t2 = int(name.split("_")[3])
            return t1, t2
        except Exception:
            return None, None

    scored = []
    for ev in casa_event_dirs:
        t1, t2 = parse_times(ev)
        if sun_start is None or t1 is None or t2 is None:
            scored.append((10**9, ev))
            continue
        # prefer before: distance = sun_start - casa_end if casa_end <= sun_start
        if t2 <= sun_start:
            dist = sun_start - t2
            scored.append((dist, ev))
        else:
            # after: penalize so it ranks after any "before"
            dist = 10**6 + (t1 - sun_start)
            scored.append((dist, ev))
    scored.sort(key=lambda x: x[0])
    return scored[0][1] if scored else None

def split_casa_candidates_relative_to_sun(sun_event_dir, casa_event_dirs):
    """
    Split CASA candidates into (before_list, after_list) relative to SUN start time.
    """
    sun_name = Path(sun_event_dir).name
    try:
        sun_start = int(sun_name.split("_")[1])  # 101000
    except Exception:
        sun_start = None

    before, after = [], []
    for ev in casa_event_dirs:
        name = Path(ev).name
        try:
            t1 = int(name.split("_")[1])
            t2 = int(name.split("_")[3])
        except Exception:
            # if parse fails, just put to 'after'
            after.append(ev)
            continue

        if sun_start is None:
            after.append(ev)
        elif t2 <= sun_start:
            before.append(ev)
        else:
            after.append(ev)

    return sorted(before), sorted(after)