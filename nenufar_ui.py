# nenufar_ui.py

from pathlib import Path
import re
import glob
import ipywidgets as w
from IPython.display import display

import datetime



def load_and_select_sb(BASE, WORKROOT):
    """
    Build interactive UI for:
      - Year / Month / Day selection
      - CAS_A_TRACKING selection
      - Frequency range selection
      - Run scan_sun_and_casa_by_ymd

    UI additions (v2 display only):
      - Show: 30/100/300/All
      - Search filter (optional)
      - Print total rows found
    """
    import re
    from pathlib import Path

    import ipywidgets as w
    from IPython.display import display

    from nenufar_sb_scan import (
        scan_sun_and_casa_by_ymd,
        list_event_dirs_by_date,
        pick_event_dir,
        pick_closest_calibrator,
    )

    BASE = Path(BASE)
    WORKROOT = str(WORKROOT)

    # ---------- build availability index ----------
    pat = re.compile(r"(\d{4})(\d{2})(\d{2})_.*_SUN_TRACKING$")

    avail = {}
    for ydir in sorted(BASE.glob("[0-9]" * 4)):
        if not ydir.is_dir():
            continue
        y = int(ydir.name)
        for mdir in sorted(ydir.glob("[0-9]" * 2)):
            if not mdir.is_dir():
                continue
            m = int(mdir.name)
            days = set()
            for ev in mdir.glob("*SUN_TRACKING"):
                mm = pat.search(ev.name)
                if mm:
                    days.add(int(mm.group(3)))
            if days:
                avail.setdefault(y, {})[m] = sorted(days)

    years = sorted(avail.keys(), reverse=True)

    def months_for(y):
        return sorted(avail.get(y, {}).keys())

    def days_for(y, m):
        return avail.get(y, {}).get(m, [])

    # ---------- widgets ----------
    year = w.Dropdown(options=years, description="Year")
    month = w.Dropdown(options=months_for(year.value), description="Month")
    day = w.Dropdown(options=days_for(year.value, month.value), description="Day")

    cal = w.Dropdown(options=[], description="Cal")

    allfreq = w.Checkbox(value=False, description="All freq")
    fmin = w.FloatText(value=70.0, description="Fmin")
    fmax = w.FloatText(value=80.0, description="Fmax")

    # display controls
    preview = w.Dropdown(options=[30, 100, 300, "All"], value=30, description="Show")
    search = w.Text(value="", description="Search", placeholder="e.g. SB407 / 407 / 60.3")

    btn = w.Button(description="Run", button_style="success")
    out = w.Output()

    # ---------- refresh calibrator ----------
    def refresh_cal():
        try:
            sun_ev = pick_event_dir(str(BASE), year.value, month.value, day.value, "SUN_TRACKING")
            cands = list_event_dirs_by_date(str(BASE), year.value, month.value, day.value, "CAS_A_TRACKING")
            cal.options = cands
            if cands:
                cal.value = pick_closest_calibrator(sun_ev, cands)
        except Exception:
            cal.options = []
            cal.value = None

    refresh_cal()

    # ---------- linkage ----------
    def on_year_change(change):
        if change.get("name") != "value":
            return
        mopts = months_for(year.value)
        month.options = mopts
        if mopts:
            month.value = mopts[0]
        refresh_cal()

    def on_month_change(change):
        if change.get("name") != "value":
            return
        dopts = days_for(year.value, month.value)
        day.options = dopts
        if dopts:
            day.value = dopts[0]
        refresh_cal()

    def on_day_change(change):
        if change.get("name") != "value":
            return
        refresh_cal()

    year.observe(on_year_change, names="value")
    month.observe(on_month_change, names="value")
    day.observe(on_day_change, names="value")

    # ---------- helpers ----------
    def _apply_search(df, q: str):
        q = (q or "").strip()
        if not q:
            return df
        qlow = q.lower()
        # simple, robust: match anywhere in row string
        mask = df.astype(str).apply(lambda row: row.str.lower().str.contains(qlow, na=False)).any(axis=1)
        return df.loc[mask]
    

    # ---------- run ----------
    def run(_):
        with out:
            out.clear_output()

            fr = None if allfreq.value else (float(fmin.value), float(fmax.value))
            df_sel, meta = scan_sun_and_casa_by_ymd(
                str(BASE),
                WORKROOT,
                year.value,
                month.value,
                day.value,
                fr,
                casa_event_dir=cal.value,
            )

            # filter
            df_show = _apply_search(df_sel, search.value)

            print(meta)
            print(f"Found {len(df_sel)} SBs total; after search filter: {len(df_show)}")
            print(f"Showing: {preview.value}")

            if preview.value == "All":
                display(df_show)
            else:
                display(df_show.head(int(preview.value)))

    btn.on_click(run)

    display(
        w.HBox([year, month, day, cal]),
        w.HBox([allfreq, fmin, fmax, preview, search, btn]),
        out,
    )

def run_step1_ui(plan_json_path, out_root=None):
    """
    Small UI to run DP3 Step-1 (prepare Cal + Sun) for selected SBs.

    Inputs
    - plan_json_path: imaging_plan.json or selected_sb_pair_list.json
        must include: selected_sb (or selected_sbs) + sun_ms
        optionally: casa_pre_ms/casa_post_ms (or casa_pre/casa_post) and/or casa_ms
    - out_root: where to write outputs; default: alongside plan json folder

    What it runs (per SB)
    - 1A: Cal  AOFlagger + averager  -> CasA_{SB}_prep.MS   (new MS)
    - 1B: Cal  preflagger (in-place) -> modifies CasA_{SB}_prep.MS
    - 1C: Sun  clear flags           -> SUN_{SB}_prep.MS    (new MS)

    New behavior
    - Overwrite checkbox:
        * OFF (default): if output already exists -> SKIP that step
        * ON: delete existing outputs first -> re-run steps
    """
    import json
    from pathlib import Path
    import ipywidgets as w
    from IPython.display import display
    import subprocess, shlex
    import shutil

    # ---- container command ----
    SIF = "/home/jzhang/LOFARimgCode/linc_latest.sif"
    BIND = "-B /databf -B /data"
    ENGINE = "apptainer"   # if you only have singularity -> "singularity"
    DP3_CMD = f"{ENGINE} exec {BIND} {SIF} DP3"

    plan_json_path = Path(plan_json_path)
    payload = json.load(open(plan_json_path, "r"))

    # ---- locate arrays ----
    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    sun_ms_list = payload.get("sun_ms") or []
    casa_ms_list = payload.get("casa_ms") or []  # might be chosen already
    casa_pre_ms_list  = payload.get("casa_pre_ms")  or payload.get("casa_pre")
    casa_post_ms_list = payload.get("casa_post_ms") or payload.get("casa_post")

    if not (sbs and sun_ms_list):
        raise ValueError("JSON must contain 'selected_sb' (or 'selected_sbs') and 'sun_ms'.")

    # build mapping sb -> sun_ms
    sun_map = dict(zip(sbs, sun_ms_list))

    # cal mode: if pre/post lists exist, allow those; else allow dropdown-only
    cal_mode_opts = ["dropdown"]
    if casa_pre_ms_list and casa_post_ms_list:
        cal_mode_opts = ["dropdown", "pre", "post"]

    cal_mode = w.ToggleButtons(options=cal_mode_opts, value="dropdown", description="Cal mode")

    sb_select = w.SelectMultiple(
        options=sbs,
        value=(("SB359.MS",) if ("SB359.MS" in sbs) else (sbs[0],)) if sbs else (),
        description="SBs",
        layout=w.Layout(width="95%", height="180px"),
    )

    overwrite = w.Checkbox(value=False, description="Overwrite existing outputs")
    run_btn = w.Button(description="Run Step-1", button_style="success")
    out = w.Output()

    if out_root is None:
        out_root = plan_json_path.parent
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    def _rm_tree(p: Path):
        if p.exists():
            shutil.rmtree(p)

    def _run_to_log(cmd: str, log_path: Path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as f:
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, text=True)
            p.wait()
            if p.returncode != 0:
                raise RuntimeError(f"Command failed (see log): {log_path}\nCMD: {cmd}")

    def get_casa_ms_for_sb(sb):
        """
        Return a CASA MS path for this SB.

        Accepts either:
        - full MS path: .../CAS_A_TRACKING/L1/SB359.MS
        - event dir:    .../CAS_A_TRACKING   (then we append /L1/SBxxx.MS)
        """
        from pathlib import Path

        def _normalize(entry):
            if entry is None:
                return None
            p = Path(str(entry))

            # case A: already an MS path (endswith .MS or exists as a directory MS)
            if str(p).endswith(".MS") or p.name.endswith(".MS"):
                return str(p)

            # case B: event dir like ..._CAS_A_TRACKING (or any dir): append L1/SBxxx.MS
            # also handle the case where entry accidentally points to .../L1
            if p.name == "L1":
                return str(p / sb)

            return str(p / "L1" / sb)

        # priority: pre/post if available and user selected, else use existing casa_ms list
        if cal_mode.value == "pre" and casa_pre_ms_list:
            entry = dict(zip(sbs, casa_pre_ms_list)).get(sb)
            return _normalize(entry)

        if cal_mode.value == "post" and casa_post_ms_list:
            entry = dict(zip(sbs, casa_post_ms_list)).get(sb)
            return _normalize(entry)

        # dropdown: use casa_ms list if present
        if casa_ms_list:
            entry = dict(zip(sbs, casa_ms_list)).get(sb)
            return _normalize(entry)

        raise ValueError("No CASA path available in JSON.")

    def on_run(_):
        with out:
            out.clear_output()

            chosen = list(sb_select.value)
            if not chosen:
                print("No SB selected.")
                return

            print("Plan:", plan_json_path)
            print("Out root:", out_root)
            print("Cal mode:", cal_mode.value)
            print("Overwrite:", overwrite.value)
            print("SBs:", ", ".join(chosen))

            for sb in chosen:
                if sb not in sun_map:
                    print(f"[WARN] SB not in sun_map: {sb} (skip)")
                    continue

                sun_ms = sun_map[sb]
                casa_ms = get_casa_ms_for_sb(sb)
                if not casa_ms:
                    print(f"[WARN] No CASA MS found for {sb} under cal_mode={cal_mode.value} (skip)")
                    continue

                sb_tag = sb.replace(".MS", "")
                sb_dir = out_root / sb_tag
                sb_dir.mkdir(parents=True, exist_ok=True)

                casa_out = sb_dir / f"CasA_{sb_tag}_prep.MS"
                sun_out  = sb_dir / f"SUN_{sb_tag}_prep.MS"

                log1 = sb_dir / "dp3_casa_aoflag_avg.log"
                log2 = sb_dir / "dp3_casa_preflag.log"
                log3 = sb_dir / "dp3_sun_clearflags.log"

                # overwrite policy: remove existing outputs first
                if overwrite.value:
                    _rm_tree(casa_out)
                    _rm_tree(sun_out)

                # 1A: Cal AOFlagger + avg (create casa_out)
                cmd1 = (
                    f"{DP3_CMD} msin={shlex.quote(casa_ms)} "
                    f"msout={shlex.quote(str(casa_out))} "
                    f"steps=[flag,averager] "
                    f"flag.type=aoflagger "
                    f"averager.timestep=1 averager.freqstep=1"
                )
                print(f"\n=== 1A: {sb} ===")
                if (not overwrite.value) and casa_out.exists():
                    print(f"[SKIP 1A] exists: {casa_out}")
                else:
                    print(cmd1)
                    _run_to_log(cmd1, log1)

                # 1B: Cal preflagger (in-place on casa_out)
                cmd2 = (
                    f"{DP3_CMD} msin={shlex.quote(str(casa_out))} "
                    f"msout=. "
                    f"steps=[flag] "
                    f"flag.type=preflagger "
                    f"flag.baseline='MR102NEN&&*;MR103NEN&&*'"
                )
                print(f"=== 1B: {sb} ===")
                if not casa_out.exists():
                    raise RuntimeError(f"1B needs {casa_out}, but it does not exist.")
                # minimal marker: if log2 exists and not overwriting -> skip
                if (not overwrite.value) and log2.exists():
                    print(f"[SKIP 1B] log exists: {log2}")
                else:
                    print(cmd2)
                    _run_to_log(cmd2, log2)

                # 1C: Sun clear flags -> create sun_out
                cmd3 = (
                    f"{DP3_CMD} msin={shlex.quote(sun_ms)} "
                    f"msout={shlex.quote(str(sun_out))} "
                    f"steps=[flag] "
                    f"flag.type=preflagger flag.mode=clear flag.baseline='*&&*'"
                )
                print(f"=== 1C: {sb} ===")
                if (not overwrite.value) and sun_out.exists():
                    print(f"[SKIP 1C] exists: {sun_out}")
                else:
                    print(cmd3)
                    _run_to_log(cmd3, log3)

                print(f"OK: {sb} -> {sb_dir}")

            print("\nAll done.")

    run_btn.on_click(on_run)

    display(w.VBox([w.HBox([cal_mode, overwrite, run_btn]), sb_select, out]))

def run_step2_zoom_ui(
    plan_json_path,
    step1_root=None,
    out_root=None,
    default_start="2024/03/10/10:10:00",
    default_end="2024/03/10/10:20:00",
):
    """
    Step-2A (Zooming in):
      Cut SUN MS to a ROI time range using DP3 msin.starttime/endtime (steps=[]).

    Inputs
    - plan_json_path: selected_sb_pair_list.json (needs selected_sb + sun_ms at least)
    - step1_root: Step-1 outputs root (e.g. .../step1_outputs_20240310). If provided and
                  SUN_{SB}_prep.MS exists, we use it as msin for ROI cutting.
    - out_root: where to write step-2 outputs. default: sibling folder of plan
    """

    import json
    from pathlib import Path
    import ipywidgets as w
    from IPython.display import display
    import subprocess, shlex

    # ---- container / dp3 ----
    SIF = "/home/jzhang/LOFARimgCode/linc_latest.sif"
    BIND = "-B /databf -B /data"
    ENGINE = "apptainer"   # if you only have singularity, change to "singularity"
    DP3_CMD = f"{ENGINE} exec {BIND} {SIF} DP3"

    plan_json_path = Path(plan_json_path)
    payload = json.load(open(plan_json_path, "r"))

    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    sun_ms_list = payload.get("sun_ms") or []
    if not (sbs and sun_ms_list):
        raise ValueError("JSON must contain 'selected_sb' (or selected_sbs) and 'sun_ms'.")

    sun_map = dict(zip(sbs, sun_ms_list))

    # defaults
    if out_root is None:
        out_root = plan_json_path.parent.parent / "step2_outputs"
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    step1_root = Path(step1_root) if step1_root else None

    # ---- widgets ----
    sb_select = w.SelectMultiple(
        options=sbs,
        value=(sbs[0],) if ("SB359.MS" in sbs) else tuple(sbs[:1]),
        description="SBs",
        layout=w.Layout(width="95%", height="180px"),
    )

    t_start = w.Text(value=default_start, description="Start", layout=w.Layout(width="420px"))
    t_end   = w.Text(value=default_end,   description="End",   layout=w.Layout(width="420px"))

    overwrite = w.Checkbox(value=False, description="Overwrite existing outputs")
    run_btn = w.Button(description="Run Step-2 (ROI)", button_style="success")
    out = w.Output()

    def _run_to_log(cmd: str, log_path: Path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as f:
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, text=True)
            p.wait()
            if p.returncode != 0:
                raise RuntimeError(f"Command failed (see log): {log_path}\nCMD: {cmd}")

    def _resolve_sun_msin(sb: str) -> str:
        """
        Prefer Step-1 prepared SUN MS:
          {step1_root}/{SBtag}/SUN_{SBtag}_prep.MS
        else fallback to plan sun_ms.
        """
        if step1_root:
            sb_tag = sb.replace(".MS", "")
            cand = step1_root / sb_tag / f"SUN_{sb_tag}_prep.MS"
            if cand.exists():
                return str(cand)
        return str(sun_map[sb])

    def _safe_remove_ms(path: Path):
        # MS is a directory; remove recursively
        import shutil
        if path.exists():
            shutil.rmtree(path)

    def on_run(_):
        with out:
            out.clear_output()

            chosen = list(sb_select.value)
            if not chosen:
                print("No SB selected.")
                return

            print("Plan:", plan_json_path)
            print("Step1 root:", str(step1_root) if step1_root else "(none)")
            print("Out root:", out_root)
            print("ROI:", t_start.value, "->", t_end.value)
            print("Overwrite:", overwrite.value)
            print("SBs:", ", ".join(chosen))

            for sb in chosen:
                sb_tag = sb.replace(".MS", "")
                sb_dir = out_root / "ROI" / sb_tag
                sb_dir.mkdir(parents=True, exist_ok=True)

                msin = _resolve_sun_msin(sb)
                msout = sb_dir / f"SUN_{sb_tag}_ROI.MS"
                log_roi = sb_dir / "dp3_roi_cut.log"

                if msout.exists() and overwrite.value:
                    _safe_remove_ms(msout)

                if msout.exists() and not overwrite.value:
                    print(f"\n[Skip] ROI exists for {sb}: {msout}")
                    continue

                cmd_roi = (
                    f"{DP3_CMD} "
                    f"msin={shlex.quote(str(msin))} "
                    f"msout={shlex.quote(str(msout))} "
                    f"steps=[] "
                    f"msin.starttime={shlex.quote(t_start.value)} "
                    f"msin.endtime={shlex.quote(t_end.value)}"
                )

                print(f"\n=== ROI cut {sb} ===")
                print(cmd_roi)
                _run_to_log(cmd_roi, log_roi)
                print("OK:", msout)

            print("\nAll done.")

    run_btn.on_click(on_run)

    ui = w.VBox([
        w.HBox([run_btn, overwrite]),
        w.HBox([t_start, t_end]),
        sb_select,
        out
    ])
    display(ui)

def run_step3_calib_ui(
    plan_json_path,
    step1_root,
    step2_root,
    out_root,
    sourcedb,
    sif="/home/jzhang/LOFARimgCode/linc_latest.sif",
    engine="apptainer",
    bind="-B /databf -B /data",
):
    """
    Step-3: Calibration (gaincal on CasA prep MS -> applycal on Sun ROI MS)

    Inputs
    - plan_json_path: selected_sb_pair_list.json (must include selected_sb / selected_sbs)
    - step1_root: Step-1 outputs root (contains SBxxx/CasA_SBxxx_prep.MS)
    - step2_root: Step-2 outputs root (contains ROI/SBxxx/SUN_SBxxx_ROI.MS)
    - out_root:  Step-3 outputs root (writes parsets + logs under SBxxx/)
    - sourcedb: existing CasA.sourcedb path (REQUIRED)

    What it does per SB
    A) gaincal on calibrator prep MS (in-place) => writes parmdb: <calib_ms>/instrument
    B) applycal on Sun ROI MS (in-place) DATA -> CORR_NO_BEAM
    """

    import json
    from pathlib import Path
    import ipywidgets as w
    from IPython.display import display
    import subprocess
    import shlex
    import time

    # ---- paths ----
    plan_json_path = Path(plan_json_path)
    step1_root = Path(step1_root)
    step2_root = Path(step2_root)
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    sourcedb = Path(sourcedb)
    if not sourcedb.exists():
        raise FileNotFoundError(f"sourcedb not found: {sourcedb}")

    DP3_CMD = f"{engine} exec {bind} {shlex.quote(sif)} DP3"

    payload = json.load(open(plan_json_path, "r"))
    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    if not sbs:
        raise ValueError("JSON must contain 'selected_sb' or 'selected_sbs'.")

    # Cal mode only chooses WHICH original CASA you intended upstream.
    # In Step-3 we ALWAYS read Step-1 products; so cal_mode is for bookkeeping consistency.
    casa_pre_list  = payload.get("casa_pre_ms")  or payload.get("casa_pre")
    casa_post_list = payload.get("casa_post_ms") or payload.get("casa_post")
    casa_ms_list   = payload.get("casa_ms")      or []

    cal_mode_opts = ["dropdown"]
    if casa_pre_list and casa_post_list:
        cal_mode_opts = ["dropdown", "pre", "post"]

    cal_mode = w.ToggleButtons(options=cal_mode_opts, value="dropdown", description="Cal mode")

    overwrite = w.Checkbox(value=False, description="Overwrite existing outputs")

    sb_select = w.SelectMultiple(
        options=sbs,
        value=(sbs[0],) if ("SB359.MS" in sbs) else tuple(sbs[:1]),
        description="SBs",
        layout=w.Layout(width="95%", height="180px"),
    )

    run_btn = w.Button(description="Run Step-3 (Calib)", button_style="success")
    out = w.Output()

    def _run_to_log(cmd: str, log_path: Path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as f:
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, text=True)
            p.wait()
            if p.returncode != 0:
                raise RuntimeError(f"Command failed (see log): {log_path}\nCMD: {cmd}")

    def _write_text(path: Path, text: str, force: bool):
        if path.exists() and (not force):
            return False
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return True

    def _get_casa_tag_for_sb(sb: str) -> str:
        # purely for printing/debug; Step-3 uses Step-1 outputs anyway
        if cal_mode.value == "pre" and casa_pre_list:
            return "pre"
        if cal_mode.value == "post" and casa_post_list:
            return "post"
        return "dropdown"

    def on_run(_):
        with out:
            out.clear_output()

            chosen = list(sb_select.value)
            if not chosen:
                print("No SB selected.")
                return

            print("Plan:", plan_json_path)
            print("Step1 root:", step1_root)
            print("Step2 root:", step2_root)
            print("Out root:", out_root)
            print("Cal mode:", cal_mode.value)
            print("Overwrite:", overwrite.value)
            print("SBs:", ", ".join(chosen))
            print("Using sourcedb:", sourcedb)
            print()

            t0 = time.time()

            for sb in chosen:
                sb_tag = sb.replace(".MS", "")
                sb_dir = out_root / sb_tag
                sb_dir.mkdir(parents=True, exist_ok=True)

                # Step-1 products
                calib_ms = step1_root / sb_tag / f"CasA_{sb_tag}_prep.MS"
                if not calib_ms.exists():
                    raise FileNotFoundError(f"Missing Step-1 calibrator MS for {sb}: {calib_ms}")

                # Step-2 products
                sun_roi_ms = step2_root / "ROI" / sb_tag / f"SUN_{sb_tag}_ROI.MS"
                if not sun_roi_ms.exists():
                    raise FileNotFoundError(f"Missing Step-2 Sun ROI MS for {sb}: {sun_roi_ms}")

                # expected parmdb
                parmdb = calib_ms / "instrument"

                # ---- parsets ----
                par_gain = sb_dir / f"{sb_tag}_gaincal.parset"
                par_apply = sb_dir / f"{sb_tag}_applycal.parset"

                gaincal_parset = f"""
msin={calib_ms}
msout=.

steps=[gaincal]

gaincal.usebeammodel=False
gaincal.solint=1
gaincal.sources=CasA
gaincal.sourcedb={sourcedb}
gaincal.onebeamperpatch=False
gaincal.caltype=diagonal
""".lstrip()

                applycal_parset = f"""
msin={sun_roi_ms}
msout=.
msin.datacolumn=DATA
msout.datacolumn=CORR_NO_BEAM

steps=[applycal]

applycal.parmdb={parmdb}
applycal.updateweights=True
""".lstrip()

                wrote_gain = _write_text(par_gain, gaincal_parset, force=overwrite.value)
                wrote_apply = _write_text(par_apply, applycal_parset, force=overwrite.value)

                if wrote_gain:
                    print(f"[WRITE parset] {par_gain}")
                if wrote_apply:
                    print(f"[WRITE parset] {par_apply}")

                # ---- run A: gaincal (optional skip if parmdb exists and not overwrite) ----
                log_gain = sb_dir / "01_gaincal.log"
                if parmdb.exists() and (not overwrite.value):
                    print(f"[SKIP gaincal] {sb}: parmdb exists ({parmdb})")
                else:
                    print(f"\n=== Step-3A gaincal: {sb} ({_get_casa_tag_for_sb(sb)}) ===")
                    cmd_gain = f"{DP3_CMD} {shlex.quote(str(par_gain))}"
                    print(cmd_gain)
                    _run_to_log(cmd_gain, log_gain)

                # ---- run B: applycal ----
                log_apply = sb_dir / "02_applycal.log"
                print(f"=== Step-3B applycal: {sb} ===")
                cmd_apply = f"{DP3_CMD} {shlex.quote(str(par_apply))}"
                print(cmd_apply)
                _run_to_log(cmd_apply, log_apply)

                print(f"OK: {sb} -> {sb_dir}\n")

            print(f"All done. Elapsed: {time.time() - t0:.1f} s")

    run_btn.on_click(on_run)
    display(w.VBox([w.HBox([cal_mode, overwrite, run_btn]), sb_select, out]))

def run_step4_wsclean_ui(
    plan_json_path,
    step2_root,
    out_root,
    *,
    sif="/home/jzhang/LOFARimgCode/linc_latest.sif",
    bind="-B /databf -B /data",
    engine="apptainer",
    default_data_column="CORR_NO_BEAM",
    default_size=(1126, 1126),
    default_scale="50asec",
    default_weight_mode="briggs",
    default_robust=0.0,
    default_mgain=0.7,
    default_mem_gb=90,
    default_niter=12000,
    default_auto_mask=3.0,
    default_auto_threshold=0.3,
    default_intervals_out=10,
    default_interval_start="",
    default_interval_end="",
):
    """
    Step-4 (Cleaning / Imaging): run WSClean on Step-2 ROI MS (SUN_<SB>_ROI.MS),
    using DATA or CORR_NO_BEAM column, with per-SB logs and skip/overwrite behavior.
    """

    import json
    from pathlib import Path
    import time
    import subprocess, shlex
    import ipywidgets as w
    from IPython.display import display

    plan_json_path = Path(plan_json_path)
    step2_root = Path(step2_root)
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    payload = json.load(open(plan_json_path, "r"))
    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    if not sbs:
        raise ValueError("JSON must contain 'selected_sb' or 'selected_sbs'.")

    # Command base (inside container)
    WSCLEAN_CMD = f"{engine} exec {bind} {sif} wsclean"

    # ---------- helpers ----------
    def _run_to_log(cmd: str, log_path: Path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as f:
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, text=True)
            p.wait()
            if p.returncode != 0:
                raise RuntimeError(f"Command failed (see log): {log_path}\nCMD: {cmd}")

    def _roi_ms_path(sb: str) -> Path:
        tag = sb.replace(".MS", "")
        # Expected step2 output layout:
        # step2_root/ROI/SBxxx/SUN_SBxxx_ROI.MS
        return step2_root / "ROI" / tag / f"SUN_{tag}_ROI.MS"

    def _image_prefix(sb: str) -> Path:
        tag = sb.replace(".MS", "")
        # Output layout:
        # out_root/SBxxx/<prefix>.* (fits, psf, dirty, model, etc.)
        return out_root / tag / tag

    def _outputs_exist(prefix: Path) -> bool:
        # WSClean typically produces files like:
        # <prefix>-image.fits, <prefix>-dirty.fits, <prefix>-psf.fits, etc.
        # Also when intervals-out > 1: <prefix>-0000-image.fits ...
        pat1 = str(prefix) + "-image.fits"
        pat2 = str(prefix) + "-0000-image.fits"
        return Path(pat1).exists() or Path(pat2).exists()

    def _parse_size(text: str):
        # accept "1126 1126" or "1126,1126" or "1126x1126"
        s = text.strip().lower().replace("x", " ").replace(",", " ")
        parts = [p for p in s.split() if p]
        if len(parts) != 2:
            raise ValueError("Size must have two integers, e.g. '1126 1126'.")
        return int(parts[0]), int(parts[1])

    def _build_cmd(
        ms_path: Path,
        name_prefix: Path,
        data_column: str,
        size_xy,
        scale_str: str,
        weight_mode: str,
        robust: float,
        mgain: float,
        mem_gb: int,
        niter: int,
        auto_mask: float,
        auto_threshold: float,
        intervals_out: int,
        interval_start: str,
        interval_end: str,
        pol: str,
        no_reorder: bool,
        no_update_model_required: bool,
    ) -> str:
        size_x, size_y = size_xy

        args = []
        args += [f"-mem {int(mem_gb)}"]
        if no_reorder:
            args += ["-no-reorder"]
        if no_update_model_required:
            args += ["-no-update-model-required"]

        args += [
            f"-mgain {float(mgain)}",
            f"-weight {weight_mode} {float(robust)}" if weight_mode == "briggs" else f"-weight {weight_mode}",
            f"-auto-mask {float(auto_mask)}",
            f"-auto-threshold {float(auto_threshold)}",
            f"-size {int(size_x)} {int(size_y)}",
            f"-scale {scale_str}",
            f"-pol {pol}",
            f"-data-column {data_column}",
            f"-intervals-out {int(intervals_out)}",
            f"-niter {int(niter)}",
            f"-name {shlex.quote(str(name_prefix))}",
            shlex.quote(str(ms_path)),
        ]

        # Optional interval selection:
        # If both provided -> -interval <start> <end>
        # (Start/end are integers in "WSClean interval index" space.)
        ist = str(interval_start).strip()
        ied = str(interval_end).strip()
        if ist and ied:
            args.insert(args.index(f"-intervals-out {int(intervals_out)}"), f"-interval {int(ist)} {int(ied)}")
        elif ist or ied:
            raise ValueError("Interval start/end must be both set or both empty.")

        return f"{WSCLEAN_CMD} " + " ".join(args)

    # ---------- UI ----------
    cal_mode = w.ToggleButtons(options=["dropdown", "pre", "post"], value="dropdown", description="Cal mode")
    # (Kept for visual consistency with previous steps; Step-4 actually uses Step-2 ROI MS only.)
    cal_mode.layout = w.Layout(width="360px")

    overwrite = w.Checkbox(value=False, description="Overwrite existing outputs")

    data_col = w.Dropdown(
        options=["CORR_NO_BEAM", "DATA"],
        value=default_data_column,
        description="Data column",
        layout=w.Layout(width="330px"),
    )

    size_box = w.Text(
        value=f"{default_size[0]} {default_size[1]}",
        description="Size (px)",
        placeholder="e.g. 1126 1126",
        layout=w.Layout(width="330px"),
    )

    scale_box = w.Text(
        value=str(default_scale),
        description="Scale",
        placeholder="e.g. 50asec or 0.56amin",
        layout=w.Layout(width="330px"),
    )

    weight_mode = w.Dropdown(
        options=["briggs", "natural", "uniform"],
        value=default_weight_mode,
        description="Weight",
        layout=w.Layout(width="330px"),
    )
    robust_box = w.FloatText(value=float(default_robust), description="Robust", layout=w.Layout(width="330px"))

    mgain_box = w.FloatText(value=float(default_mgain), description="mgain", layout=w.Layout(width="330px"))
    mem_box = w.IntText(value=int(default_mem_gb), description="mem (GB)", layout=w.Layout(width="330px"))
    niter_box = w.IntText(value=int(default_niter), description="niter", layout=w.Layout(width="330px"))
    amask_box = w.FloatText(value=float(default_auto_mask), description="auto-mask", layout=w.Layout(width="330px"))
    athr_box = w.FloatText(value=float(default_auto_threshold), description="auto-threshold", layout=w.Layout(width="330px"))

    intervals_out_box = w.IntText(
        value=int(default_intervals_out), description="intervals-out", layout=w.Layout(width="330px")
    )
    istart_box = w.Text(value=str(default_interval_start), description="interval start", layout=w.Layout(width="330px"))
    iend_box = w.Text(value=str(default_interval_end), description="interval end", layout=w.Layout(width="330px"))

    pol_box = w.Dropdown(options=["I"], value="I", description="Pol", layout=w.Layout(width="330px"))

    no_reorder = w.Checkbox(value=True, description="-no-reorder")
    no_update_model_required = w.Checkbox(value=True, description="-no-update-model-required")

    sb_select = w.SelectMultiple(
        options=sbs,
        value=(sbs[0],),
        description="SBs",
        layout=w.Layout(width="95%", height="220px"),
    )

    run_btn = w.Button(description="Run Step-4 (WSClean)", button_style="success")
    out = w.Output()

    def _toggle_robust(*_):
        robust_box.disabled = (weight_mode.value != "briggs")

    weight_mode.observe(_toggle_robust, names="value")
    _toggle_robust()

    def on_run(_):
        t0 = time.time()
        with out:
            out.clear_output()

            chosen = list(sb_select.value)
            if not chosen:
                print("No SB selected.")
                return

            print("Plan:", plan_json_path)
            print("Step2 root:", step2_root)
            print("Out root:", out_root)
            print("Data column:", data_col.value)
            print("Overwrite:", overwrite.value)
            print("SBs:", ", ".join(chosen))
            print()

            # Parse size once
            try:
                size_xy = _parse_size(size_box.value)
            except Exception as e:
                raise ValueError(f"Invalid size: {e}")

            for sb in chosen:
                tag = sb.replace(".MS", "")
                ms_in = _roi_ms_path(sb)
                prefix = _image_prefix(sb)
                sb_dir = prefix.parent
                sb_dir.mkdir(parents=True, exist_ok=True)

                if not ms_in.exists():
                    print(f"[SKIP] {sb}: ROI MS not found: {ms_in}")
                    continue

                if (not overwrite.value) and _outputs_exist(prefix):
                    print(f"[SKIP] {sb}: outputs exist (enable overwrite to rerun): {prefix}*")
                    continue

                log_path = sb_dir / "wsclean.log"

                cmd = _build_cmd(
                    ms_path=ms_in,
                    name_prefix=prefix,
                    data_column=data_col.value,
                    size_xy=size_xy,
                    scale_str=scale_box.value.strip(),
                    weight_mode=weight_mode.value,
                    robust=robust_box.value,
                    mgain=mgain_box.value,
                    mem_gb=mem_box.value,
                    niter=niter_box.value,
                    auto_mask=amask_box.value,
                    auto_threshold=athr_box.value,
                    intervals_out=intervals_out_box.value,
                    interval_start=istart_box.value,
                    interval_end=iend_box.value,
                    pol=pol_box.value,
                    no_reorder=no_reorder.value,
                    no_update_model_required=no_update_model_required.value,
                )

                print(f"=== Step-4 WSClean: {sb} ===")
                print(cmd)
                _run_to_log(cmd, log_path)
                print(f"OK: {sb} -> {sb_dir}")
                print(f"Log: {log_path}")
                print()

            print(f"All done. Elapsed: {time.time() - t0:.1f} s")

    run_btn.on_click(on_run)

    # Layout
    left = w.VBox([cal_mode, sb_select])
    mid = w.VBox(
        [
            data_col,
            size_box,
            scale_box,
            weight_mode,
            robust_box,
            mgain_box,
            mem_box,
            niter_box,
            amask_box,
            athr_box,
        ],
        layout=w.Layout(width="360px"),
    )
    right = w.VBox(
        [
            intervals_out_box,
            istart_box,
            iend_box,
            pol_box,
            w.HBox([no_reorder, no_update_model_required]),
        ],
        layout=w.Layout(width="360px"),
    )

    header = w.HBox([overwrite, run_btn])
    display(w.VBox([header, w.HBox([left, mid, right]), out]))

def run_step4_quicklook_ui(
    step4_root,
    out_root=None,
    default_sb="SB359",
    default_crop_half_width_arcsec=1000,
    default_clim_pct=(5, 99),
    default_contours=(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95),
    site_lat_deg=47.382,
    site_lon_deg=2.195,
    site_height_m=136.0,
):
    """
    Step-4 Quicklook UI:
      - scan step4_root/SBxxx/*-image.fits
      - multi-select FITS
      - generate png quicklooks (+ optional video)

    Input:
      step4_root: e.g. /data/jzhang/nenufar_workflows/step4_outputs_20240310
      out_root:   where to save quicklooks; default step4_root/quicklook

    Output:
      out_root/SBxxx/*.png
      out_root/SBxxx/movie.mp4 (optional)
      out_root/SBxxx/quicklook.log
    """
    import os
    import re
    import glob
    import shutil
    import subprocess
    import numpy as np
    import ipywidgets as w
    from pathlib import Path
    from IPython.display import display

    # ---------- helpers ----------
    def _sb_tag(sb):
        return sb.replace(".MS", "").replace("SB", "SB")

    def _find_sbs(step4_root_path: Path):
        sbs = []
        for p in sorted(step4_root_path.glob("SB*")):
            if p.is_dir() and re.match(r"SB\d{3}$", p.name):
                sbs.append(p.name + ".MS")
        return sbs

    def _list_fits(step4_root_path: Path, sb: str):
        tag = _sb_tag(sb)  # "SB359"
        sb_dir = step4_root_path / tag
        pats = [
            str(sb_dir / f"{tag}-t*-image.fits"),
            str(sb_dir / f"{tag}_t*-image.fits"),
            str(sb_dir / "*-image.fits"),
        ]
        files = []
        for pat in pats:
            files.extend(glob.glob(pat))
        files = sorted(set(files))
        return files

    def _parse_tindex(fname: str):
        # SB359-t0007-image.fits -> 7
        m = re.search(r"-t(\d+)-", os.path.basename(fname))
        return int(m.group(1)) if m else 10**9

    def _ensure_dir(p: Path):
        p.mkdir(parents=True, exist_ok=True)
        return p

    def _run_to_log(cmd, log_path: Path):
        log_path = Path(log_path)
        with open(log_path, "a") as f:
            f.write("\n" + "=" * 80 + "\n")
            f.write(cmd + "\n")
            f.flush()
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, text=True)
            ret = p.wait()
        if ret != 0:
            raise RuntimeError(f"Command failed (see log): {log_path}\nCMD: {cmd}")

    def _load_sort_fits(files):
        # sort by t-index then name
        return sorted(files, key=lambda x: (_parse_tindex(x), os.path.basename(x)))
    

    # ---------- plotting core ----------
    def _quicklook_one(
        fits_path: str,
        out_png: str,
        crop_half_width_arcsec: float,
        clim_pct=(5, 99),
        contours_perc=None,
        site=None,
        overwrite=False,
        draw_beam=True,
        cmap="viridis",
    ):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter

        import astropy.units as u
        from astropy.io import fits
        from astropy.time import Time
        from astropy.coordinates import SkyCoord

        import sunpy.map
        from sunpy.coordinates import frames, sun

        out_png = str(out_png)
        if (not overwrite) and os.path.exists(out_png):
            return out_png

        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data

        # ensure 2D
        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image: {fits_path} shape={getattr(data, 'shape', None)}")

        # metadata
        obstime = None
        for k in ["DATE-OBS", "DATEOBS", "DATE_OBS", "DATE"]:
            if k in hdr:
                try:
                    obstime = Time(hdr[k])
                    break
                except Exception:
                    pass
        if obstime is None:
            # fallback: now
            obstime = Time.now()

        freq_Hz = hdr.get("CRVAL3", None)
        frequency = (freq_Hz * u.Hz) if freq_Hz is not None else (np.nan * u.Hz)

        # pixel scale (deg/pix -> arcsec/pix), positive
        cdelt1 = abs(hdr.get("CDELT1", np.nan)) * u.deg
        cdelt2 = abs(hdr.get("CDELT2", np.nan)) * u.deg
        cdelt1 = cdelt1.to(u.arcsec) if np.isfinite(cdelt1.value) else (np.nan * u.arcsec)
        cdelt2 = cdelt2.to(u.arcsec) if np.isfinite(cdelt2.value) else (np.nan * u.arcsec)

        # ---- try SunPy Map in Helioprojective (best) ----
        rsub = None
        used_sunpy = True
        try:
            # observer at site (GCRS)
            site_gcrs = SkyCoord(site.get_gcrs(obstime))

            # reference sky coord from header CRVAL1/2
            import astropy.units as u
            cunit1 = u.Unit(hdr.get("CUNIT1", "deg"))
            cunit2 = u.Unit(hdr.get("CUNIT2", "deg"))

            ref_gcrs = SkyCoord(
                hdr["CRVAL1"] * cunit1,
                hdr["CRVAL2"] * cunit2,
                frame="gcrs",
                obstime=obstime,
                obsgeoloc=site_gcrs.cartesian,
                obsgeovel=site_gcrs.velocity.to_cartesian(),
                distance=site_gcrs.hcrs.distance,
            )

            ref_hpc = ref_gcrs.transform_to(frames.Helioprojective(observer=site_gcrs))

            # rotate so solar north up
            P1 = sun.P(obstime)

            # reference pixel is 0-based in make_fitswcs_header
            ref_pix = np.array([hdr["CRPIX1"] - 1, hdr["CRPIX2"] - 1]) * u.pixel
            scale = np.array([cdelt1.value, cdelt2.value]) * (u.arcsec / u.pixel)

            new_header = sunpy.map.make_fitswcs_header(
                data=data,
                coordinate=ref_hpc,
                reference_pixel=ref_pix,
                scale=scale,
                rotation_angle=-P1,
                wavelength=frequency.to(u.MHz) if np.isfinite(frequency.value) else None,
                observatory="NenuFAR (Nançay)",
            )

            rmap = sunpy.map.Map(data, new_header)
            rmap_rot = rmap.rotate()

            # crop
            hw = float(crop_half_width_arcsec) * u.arcsec
            bl = SkyCoord(-hw, -hw, frame=rmap_rot.coordinate_frame)
            tr = SkyCoord(hw, hw, frame=rmap_rot.coordinate_frame)
            rsub = rmap_rot.submap(bl, top_right=tr)
        except Exception:
            used_sunpy = False
            rsub = None

        # ---- plot layout ----
        fig = plt.figure(figsize=(8.8, 7.4), constrained_layout=False)
        gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[1.0, 0.04], wspace=0.15)

        if used_sunpy and rsub is not None:
            ax = fig.add_subplot(gs[0, 0], projection=rsub)
            cax = fig.add_subplot(gs[0, 1])

            vmin, vmax = np.nanpercentile(rsub.data, clim_pct)
            im = rsub.plot(axes=ax, cmap=cmap, vmin=vmin, vmax=vmax)

            try:
                rsub.draw_limb(axes=ax)
            except Exception:
                pass
            try:
                rsub.draw_grid(axes=ax)
            except Exception:
                pass

            if contours_perc:
                try:
                    import astropy.units as u
                    rsub.draw_contours(np.array(contours_perc) * u.percent, colors="k",
                                       linewidths=1.1, alpha=0.8, axes=ax)
                except Exception:
                    pass

        else:
            # fallback: pixel plot
            ax = fig.add_subplot(gs[0, 0])
            cax = fig.add_subplot(gs[0, 1])

            vmin, vmax = np.nanpercentile(data, clim_pct)
            im = ax.imshow(data, origin="lower", vmin=vmin, vmax=vmax)
            ax.set_xlabel("X (pix)")
            ax.set_ylabel("Y (pix)")

        # colorbar
        cbar = fig.colorbar(im, cax=cax)
        sf = ScalarFormatter(useMathText=True)
        sf.set_powerlimits((-2, 3))
        cbar.formatter = sf
        cbar.update_ticks()
        cbar.set_label("Stokes I (arb.)", labelpad=12)

        # beam (optional)
        if draw_beam:
            try:
                from matplotlib.patches import Ellipse
                bmaj = hdr.get("BMAJ"); bmin = hdr.get("BMIN"); bpa = hdr.get("BPA")
                if bmaj is not None and bmin is not None and bpa is not None:
                    # in arcsec
                    import astropy.units as u
                    bmaj_as = (abs(bmaj) * u.deg).to_value(u.arcsec)
                    bmin_as = (abs(bmin) * u.deg).to_value(u.arcsec)
                    # place in lower-left (data coords)
                    x0 = ax.get_xlim()[0] + 0.08*(ax.get_xlim()[1]-ax.get_xlim()[0])
                    y0 = ax.get_ylim()[0] + 0.08*(ax.get_ylim()[1]-ax.get_ylim()[0])
                    e = Ellipse((x0, y0), width=bmaj_as, height=bmin_as, angle=-float(bpa),
                                edgecolor="w", facecolor="none", lw=1.5)
                    ax.add_patch(e)
            except Exception:
                pass

        # titles
        left_title = f"{os.path.basename(fits_path)}   {obstime.isot}"
        fig.text(0.02, 0.98, left_title, ha="left", va="top", fontsize=16)

        if np.isfinite(frequency.value):
            freq_str = f"{frequency.to_value(u.MHz):.1f} MHz"
            fig.text(0.98, 0.98, freq_str, ha="right", va="top",
                     fontsize=15, bbox=dict(fc=(1, 1, 1, 0.7), ec="0.7", pad=2))

        _ensure_dir(Path(out_png).parent)
        plt.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_png

    def _make_video_ffmpeg(png_list, out_mp4, fps=8):
        """
        Robust mp4 maker:
        - avoid concat demuxer timestamp issues (DTS/PTS invalid dropping)
        - ensure even width/height for libx264 (e.g., 1298x1061 -> 1298x1060)
        """
        from pathlib import Path
        import os
        import shutil
        import subprocess

        png_list = list(png_list)
        if len(png_list) == 0:
            return None

        out_mp4 = Path(out_mp4)

        # ---- build a numbered image sequence via symlinks ----
        # e.g. _frames/SB359_000000.png, SB359_000001.png, ...
        frames_dir = out_mp4.parent / (out_mp4.stem + "_frames")
        if frames_dir.exists():
            shutil.rmtree(frames_dir)
        frames_dir.mkdir(parents=True, exist_ok=True)

        # keep deterministic order
        png_list = sorted(png_list)

        for i, p in enumerate(png_list):
            p = Path(p)
            link_path = frames_dir / f"{i:06d}.png"
            try:
                os.symlink(p, link_path)   # fast
            except Exception:
                shutil.copy2(p, link_path)  # fallback if symlink not allowed

        # ---- ffmpeg: image sequence -> mp4 ----
        # scale=trunc(iw/2)*2:trunc(ih/2)*2 ensures even dims for x264
        cmd = [
            "ffmpeg",
            "-y",
            "-framerate", str(int(fps)),
            "-i", str(frames_dir / "%06d.png"),
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-movflags", "+faststart",
            str(out_mp4),
        ]

        subprocess.check_call(cmd)

        # optional: cleanup frames dir (comment out if you want to keep it for debugging)
        # shutil.rmtree(frames_dir, ignore_errors=True)

        return str(out_mp4)

    def _make_video_imageio(png_list, out_gif, fps=8):
        import imageio.v2 as imageio
        imgs = []
        for p in png_list:
            imgs.append(imageio.imread(p))
        out_gif = str(out_gif)
        imageio.mimsave(out_gif, imgs, duration=1.0/max(fps, 1))
        return out_gif

    # ---------- UI ----------
    step4_root = Path(step4_root)
    if out_root is None:
        out_root = step4_root / "quicklook"
    out_root = Path(out_root)

    sbs = _find_sbs(step4_root)
    if len(sbs) == 0:
        raise ValueError(f"No SBxxx directories found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    sb_dd = w.Dropdown(options=sbs, value=default_sb_ms, description="SB")
    refresh_btn = w.Button(description="Refresh FITS list", button_style="")
    fits_sel = w.SelectMultiple(options=[], description="FITS", layout=w.Layout(width="95%", height="220px"))

    overwrite_cb = w.Checkbox(value=False, description="Overwrite PNG if exists")
    make_video_cb = w.Checkbox(value=False, description="Make video (mp4/gif)")
    fps_in = w.IntText(value=8, description="FPS")

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop half-width (arcsec)")
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low %")
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high %")
    contours_in = w.Text(value=",".join(map(str, default_contours)), description="Contours %")
    draw_beam_cb = w.Checkbox(value=True, description="Draw beam")
    run_btn = w.Button(description="Generate Quicklooks", button_style="success")

    out = w.Output()

    # site
    from astropy.coordinates import EarthLocation
    import astropy.units as u
    site = EarthLocation(lat=float(site_lat_deg)*u.deg, lon=float(site_lon_deg)*u.deg, height=float(site_height_m)*u.m)

    def _refresh(_=None):
        files = _list_fits(step4_root, sb_dd.value)
        files = _load_sort_fits(files)
        fits_sel.options = files
        with out:
            out.clear_output()
            print(f"Step4 root: {step4_root}")
            print(f"SB: {sb_dd.value}  -> found {len(files)} image.fits")
            if len(files) > 0:
                print("First:", os.path.basename(files[0]))
                print("Last :", os.path.basename(files[-1]))

    def _parse_contours(s: str):
        s = (s or "").strip()
        if not s:
            return []
        vals = []
        for x in s.split(","):
            x = x.strip()
            if not x:
                continue
            try:
                vals.append(float(x))
            except Exception:
                pass
        return vals

    def _on_run(_):
        with out:
            out.clear_output()

            sb = sb_dd.value
            files = list(fits_sel.value)
            if len(files) == 0:
                print("No FITS selected. Select one or more in the list.")
                return

            tag = _sb_tag(sb)
            sb_out = _ensure_dir(out_root / tag)
            log_path = sb_out / "quicklook.log"
            print("SB:", sb)
            print("Selected:", len(files))
            print("Out:", sb_out)
            print("Overwrite:", overwrite_cb.value)
            print("Make video:", make_video_cb.value, "FPS:", fps_in.value)

            # parse params
            clim = (float(clim_low.value), float(clim_high.value))
            contours = _parse_contours(contours_in.value)
            crop_hw = float(crop_in.value)

            pngs = []
            for fpath in files:
                out_png = sb_out / (Path(fpath).name.replace("-image.fits", "_quicklook.png").replace(".fits", "_quicklook.png"))
                try:
                    png = _quicklook_one(
                        fits_path=fpath,
                        out_png=str(out_png),
                        crop_half_width_arcsec=crop_hw,
                        clim_pct=clim,
                        contours_perc=contours,
                        site=site,
                        overwrite=overwrite_cb.value,
                        draw_beam=draw_beam_cb.value,
                    )
                    pngs.append(png)
                    with open(log_path, "a") as f:
                        f.write(f"OK  {fpath} -> {png}\n")
                except Exception as e:
                    with open(log_path, "a") as f:
                        f.write(f"FAIL {fpath} : {e}\n")
                    print("[FAIL]", os.path.basename(fpath), ":", e)

            print(f"PNG done: {len(pngs)}/{len(files)}")
            print("Log:", log_path)

            if make_video_cb.value and len(pngs) > 1:
                # try ffmpeg first
                mp4 = sb_out / f"{tag}_quicklook.mp4"
                gif = sb_out / f"{tag}_quicklook.gif"
                try:
                    shutil.which("ffmpeg")  # None if not found
                    if shutil.which("ffmpeg"):
                        _make_video_ffmpeg(pngs, mp4, fps=int(fps_in.value))
                        print("Video:", mp4)
                    else:
                        raise RuntimeError("ffmpeg not found")
                except Exception:
                    # fallback gif
                    try:
                        _make_video_imageio(pngs, gif, fps=int(fps_in.value))
                        print("GIF:", gif)
                    except Exception as e:
                        print("[WARN] video failed:", e)

    refresh_btn.on_click(_refresh)
    run_btn.on_click(_on_run)

    _refresh()

    display(
        w.VBox([
            w.HBox([sb_dd, refresh_btn, overwrite_cb, make_video_cb, fps_in]),
            fits_sel,
            w.HBox([crop_in, clim_low, clim_high]),
            w.HBox([contours_in, draw_beam_cb]),
            run_btn,
            out
        ])
    )

def run_step5a_iocorrect_solve_ui(
    step4_root,
    out_root=None,
    default_sb="SB359",
    default_crop_half_width_arcsec=5000.0,
    default_clim_pct=(5, 99),
    default_roi_arcsec=(-1500.0, 1000.0, 0.0, 1500.0),  # xmin,xmax,ymin,ymax
    thresh_frac=0.5,
    min_points=30,
    site_lat_deg=47.382,
    site_lon_deg=2.195,
    site_height_m=136.0,
):
    """
    Step-5A (IO-correct solve) UI:
      - scan step4_root/SBxxx/*-image.fits
      - choose ONE quiet-sun FITS + ROI (arcsec in rotated HPC)
      - compute centroid (2D Gaussian on ROI)
      - convert centroid -> dx_pix, dy_pix on the CROPPED map grid (centroid -> (0,0))
      - save solution JSON (for Step-5B)

    Output:
      out_root/SBxxx/step5a_solution.json
      out_root/SBxxx/step5a_quiet_preview.png  (optional: also saved)
    """
    import os
    import re
    import json
    import glob
    import datetime
    import numpy as np
    import ipywidgets as w
    from pathlib import Path
    from IPython.display import display, clear_output

    import matplotlib.pyplot as plt

    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.io import fits
    from astropy.time import Time

    import sunpy.map
    from sunpy.coordinates import frames, sun

    # scipy only used for gaussian fit
    from scipy.optimize import curve_fit



    # ------------------------
    # paths
    # ------------------------
 
    def _infer_event_tag_from_step4_root(step4_root: str | Path) -> str:
        """
        Infer event date tag like '20240310' from a step4 root path such as:
        .../step4_outputs_20240310
        Fallback to today's date if not found.
        """
        s = str(step4_root)
        m = re.search(r"(20\d{6})", s)
        if m:
            return m.group(1)
        return datetime.datetime.now().strftime("%Y%m%d")


    def _default_out_root(workroot: str | Path, prefix: str, event_tag: str) -> Path:
        """
        Standardized out_root creator:
        /data/jzhang/nenufar_workflows/{prefix}_{event_tag}
        """
        return Path(workroot) / f"{prefix}_{event_tag}"
    
    step4_root = Path(step4_root).expanduser()

    if out_root is None:
        event_tag = _infer_event_tag_from_step4_root(step4_root)
        out_root = _default_out_root("/data/jzhang/nenufar_workflows", "step_iocorrect_outputs", event_tag)

    out_root = Path(out_root).expanduser()


    # ------------------------
    # site
    # ------------------------
    site = EarthLocation(
        lat=float(site_lat_deg) * u.deg,
        lon=float(site_lon_deg) * u.deg,
        height=float(site_height_m) * u.m,
    )

    # ------------------------
    # helpers: SB + FITS list (same spirit as step4)
    # ------------------------
    def _sb_tag(sb):
        return sb.replace(".MS", "").replace("SB", "SB")

    def _find_sbs(root: Path):
        sbs = []
        for p in sorted(root.glob("SB*")):
            if p.is_dir() and re.match(r"SB\d{3}$", p.name):
                sbs.append(p.name + ".MS")
        return sbs

    def _list_fits(root: Path, sb: str):
        tag = _sb_tag(sb)
        sb_dir = root / tag
        pats = [
            str(sb_dir / f"{tag}-t*-image.fits"),
            str(sb_dir / f"{tag}_t*-image.fits"),
            str(sb_dir / "*-image.fits"),
        ]
        files = []
        for pat in pats:
            files.extend(glob.glob(pat))
        files = sorted(set(files))
        return files

    def _ensure_dir(p: Path):
        p.mkdir(parents=True, exist_ok=True)
        return p

    # ------------------------
    # core: WSClean FITS -> rotated HPC map (match Step4)
    # ------------------------
    def _read_2d_data_and_header(fpath: str):
        with fits.open(fpath) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image: {fpath} shape={getattr(data,'shape',None)}")
        return np.array(data, dtype=float), hdr

    def _build_rot_hpc_map(fpath: str):
        data, hdr = _read_2d_data_and_header(fpath)

        # obstime
        obstime = None
        for k in ["DATE-OBS", "DATEOBS", "DATE_OBS", "DATE"]:
            if k in hdr:
                try:
                    obstime = Time(hdr[k])
                    break
                except Exception:
                    pass
        if obstime is None:
            obstime = Time.now()

        # freq
        freq_Hz = hdr.get("CRVAL3", None)
        frequency = (freq_Hz * u.Hz) if freq_Hz is not None else (np.nan * u.Hz)

        # pixel scale
        cdelt1 = abs(hdr.get("CDELT1", np.nan)) * u.deg
        cdelt2 = abs(hdr.get("CDELT2", np.nan)) * u.deg
        cdelt1 = cdelt1.to(u.arcsec) if np.isfinite(cdelt1.value) else (np.nan * u.arcsec)
        cdelt2 = cdelt2.to(u.arcsec) if np.isfinite(cdelt2.value) else (np.nan * u.arcsec)

        site_gcrs = SkyCoord(site.get_gcrs(obstime))

        cunit1 = u.Unit(hdr.get("CUNIT1", "deg"))
        cunit2 = u.Unit(hdr.get("CUNIT2", "deg"))
        ref_gcrs = SkyCoord(
            hdr["CRVAL1"] * cunit1,
            hdr["CRVAL2"] * cunit2,
            frame="gcrs",
            obstime=obstime,
            obsgeoloc=site_gcrs.cartesian,
            obsgeovel=site_gcrs.velocity.to_cartesian(),
            distance=site_gcrs.hcrs.distance,
        )
        ref_hpc = ref_gcrs.transform_to(frames.Helioprojective(observer=site_gcrs))

        P1 = sun.P(obstime)

        ref_pix = np.array([hdr["CRPIX1"] - 1, hdr["CRPIX2"] - 1]) * u.pixel
        scale = np.array([cdelt1.value, cdelt2.value]) * (u.arcsec / u.pixel)

        new_header = sunpy.map.make_fitswcs_header(
            data=data,
            coordinate=ref_hpc,
            reference_pixel=ref_pix,
            scale=scale,
            rotation_angle=-P1,  # IMPORTANT: step4 style
            wavelength=frequency.to(u.MHz) if np.isfinite(frequency.value) else None,
            observatory="NenuFAR (Nançay)",
        )

        rmap = sunpy.map.Map(data, new_header)
        rmap_rot = rmap.rotate()  # IMPORTANT
        return rmap_rot, hdr, obstime, frequency

    def _submap_centered(m_rot, crop_hw_as: float):
        hw = float(crop_hw_as) * u.arcsec
        bl = SkyCoord(-hw, -hw, frame=m_rot.coordinate_frame)
        tr = SkyCoord(hw, hw, frame=m_rot.coordinate_frame)
        try:
            return m_rot.submap(bl, top_right=tr)
        except Exception:
            return m_rot

    def _submap_from_roi(m_sub, roi):
        p0 = SkyCoord(float(roi["xmin"]) * u.arcsec, float(roi["ymin"]) * u.arcsec, frame=m_sub.coordinate_frame)
        p1 = SkyCoord(float(roi["xmax"]) * u.arcsec, float(roi["ymax"]) * u.arcsec, frame=m_sub.coordinate_frame)
        return m_sub.submap(p0, top_right=p1)

    # ------------------------
    # gaussian centroid in ROI
    # ------------------------
    def _centroid_gauss(roi_map, thresh_frac=0.5, min_points=30):
        z2 = np.array(roi_map.data, dtype=float)
        if np.isnan(z2).any():
            finite = np.isfinite(z2)
            z2[~finite] = np.nanmin(z2[finite]) if finite.any() else 0.0

        amp0 = np.nanmax(z2)
        if not np.isfinite(amp0) or amp0 <= 0:
            return None

        thr = float(thresh_frac) * float(amp0)
        mask = z2 > thr
        if np.count_nonzero(mask) < int(min_points):
            return None

        yy, xx = np.indices(z2.shape)
        x = xx[mask].ravel().astype(float)
        y = yy[mask].ravel().astype(float)
        z = z2[mask].ravel().astype(float)

        def gauss2d(coords, A, x0, y0, theta, sx, sy):
            x_, y_ = coords
            ct, st = np.cos(theta), np.sin(theta)
            xp = ct * (x_ - x0) + st * (y_ - y0)
            yp = -st * (x_ - x0) + ct * (y_ - y0)
            return (A * np.exp(-0.5 * ((xp / sx) ** 2 + (yp / sy) ** 2))).ravel()

        y0_idx, x0_idx = np.unravel_index(np.nanargmax(z2), z2.shape)
        x0g, y0g = float(x0_idx), float(y0_idx)
        s_guess = max(2.0, 0.15 * min(z2.shape))
        p0 = [amp0, x0g, y0g, 0.0, s_guess, s_guess]
        bounds = (
            [0.0, 0.0, 0.0, -np.pi, 1.0, 1.0],
            [3 * amp0, z2.shape[1] - 1, z2.shape[0] - 1, np.pi, max(z2.shape), max(z2.shape)],
        )

        try:
            popt, pcov = curve_fit(gauss2d, (x, y), z, p0=p0, bounds=bounds, maxfev=20000)
        except Exception:
            return None

        A, x0_fit, y0_fit, theta_fit, sx_fit, sy_fit = popt
        cen_world = roi_map.pixel_to_world(x0_fit * u.pix, y0_fit * u.pix)

        try:
            Tx_as = float(cen_world.Tx.to_value(u.arcsec))
            Ty_as = float(cen_world.Ty.to_value(u.arcsec))
        except Exception:
            Tx_as = float(cen_world.spherical.lon.to_value(u.arcsec))
            Ty_as = float(cen_world.spherical.lat.to_value(u.arcsec))

        return dict(
            cen_world=cen_world,
            cen_Tx_as=Tx_as,
            cen_Ty_as=Ty_as,
            popt=[float(v) for v in popt],
        )

    # ------------------------
    # plot preview (same as your Step-5A quicklook you showed)
    # ------------------------
    def _plot_quiet_preview(fpath, roi, crop_hw_as, clim_pct):
        m_rot, hdr, obstime, freq = _build_rot_hpc_map(fpath)
        m_sub = _submap_centered(m_rot, crop_hw_as)
        roi_map = _submap_from_roi(m_sub, roi)

        fit = _centroid_gauss(roi_map, thresh_frac=thresh_frac, min_points=min_points)

        vmin, vmax = np.nanpercentile(m_sub.data, clim_pct)

        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(1, 2, 1, projection=m_sub)
        ax2 = fig.add_subplot(1, 2, 2, projection=roi_map)

        # left
        m_sub.plot(axes=ax1, cmap="viridis", vmin=vmin, vmax=vmax)
        try:
            m_sub.draw_limb(axes=ax1); m_sub.draw_grid(axes=ax1)
        except Exception:
            pass
        ax1.set_title(f"Full (cropped) + ROI\n{Path(fpath).name}")

        # ROI box on left
        p0 = SkyCoord(roi["xmin"] * u.arcsec, roi["ymin"] * u.arcsec, frame=m_sub.coordinate_frame)
        p1 = SkyCoord(roi["xmax"] * u.arcsec, roi["ymax"] * u.arcsec, frame=m_sub.coordinate_frame)
        x0, y0 = m_sub.world_to_pixel(p0)
        x1, y1 = m_sub.world_to_pixel(p1)
        x0 = float(np.atleast_1d(getattr(x0, "value", x0))[0])
        y0 = float(np.atleast_1d(getattr(y0, "value", y0))[0])
        x1 = float(np.atleast_1d(getattr(x1, "value", x1))[0])
        y1 = float(np.atleast_1d(getattr(y1, "value", y1))[0])
        ax1.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], transform=ax1.get_transform("pixel"), lw=2)

        # right
        roi_map.plot(axes=ax2, cmap="viridis", vmin=vmin, vmax=vmax)
        try:
            roi_map.draw_limb(axes=ax2); roi_map.draw_grid(axes=ax2)
        except Exception:
            pass
        ax2.set_title("ROI + markers")

        # disk center marker
        c0 = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=roi_map.coordinate_frame)
        ax2.plot_coord(c0, marker="x", markersize=10, mew=2)

        if fit is not None:
            ax2.plot_coord(fit["cen_world"], marker="+", markersize=14, mew=2)
            ax2.text(
                0.02, 0.98,
                f'centroid = ({fit["cen_Tx_as"]:.1f}", {fit["cen_Ty_as"]:.1f}")',
                transform=ax2.transAxes, ha="left", va="top",
                bbox=dict(fc=(1, 1, 1, 0.6), ec="none", pad=2),
            )

        plt.tight_layout()
        plt.show()

        return m_sub, fit, obstime, freq

    # ------------------------
    # compute dx_pix, dy_pix on cropped grid
    # ------------------------
    def _solve_shift(m_sub, fit):
        if fit is None:
            raise RuntimeError("Gaussian centroid fit failed. Check ROI or thresholds.")

        cen_world = fit["cen_world"]
        px_c, py_c = m_sub.world_to_pixel(cen_world)
        px_c = float(np.atleast_1d(getattr(px_c, "value", px_c))[0])
        py_c = float(np.atleast_1d(getattr(py_c, "value", py_c))[0])

        target_world = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=m_sub.coordinate_frame)
        px_t, py_t = m_sub.world_to_pixel(target_world)
        px_t = float(np.atleast_1d(getattr(px_t, "value", px_t))[0])
        py_t = float(np.atleast_1d(getattr(py_t, "value", py_t))[0])

        dx_pix = px_t - px_c
        dy_pix = py_t - py_c
        return float(dx_pix), float(dy_pix), (px_c, py_c), (px_t, py_t)

    # ------------------------
    # UI widgets
    # ------------------------
    sbs = _find_sbs(step4_root)
    if len(sbs) == 0:
        raise ValueError(f"No SBxxx directories found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    sb_dd = w.Dropdown(options=sbs, value=default_sb_ms, description="SB")
    refresh_btn = w.Button(description="Refresh FITS", button_style="")
    quiet_dd = w.Dropdown(options=[], description="Quiet FITS", layout=w.Layout(width="800px"))

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop hw (as)")
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low%")
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high%")

    xmin0, xmax0, ymin0, ymax0 = default_roi_arcsec
    xmin_w = w.FloatText(value=float(xmin0), description="xmin")
    xmax_w = w.FloatText(value=float(xmax0), description="xmax")
    ymin_w = w.FloatText(value=float(ymin0), description="ymin")
    ymax_w = w.FloatText(value=float(ymax0), description="ymax")

    run_btn = w.Button(description="Solve Step-5A (quiet centroid)", button_style="warning")
    out = w.Output()

    def _refresh(_=None):
        files = _list_fits(step4_root, sb_dd.value)
        quiet_dd.options = files
        if len(files) > 0:
            quiet_dd.value = files[0]
        with out:
            clear_output()
            print("Step4 root:", step4_root)
            print("Out root  :", out_root)
            print("SB:", sb_dd.value, "->", len(files), "image.fits")
            if len(files) > 0:
                print("First:", os.path.basename(files[0]))
                print("Last :", os.path.basename(files[-1]))

    def _on_run(_=None):
        with out:
            clear_output()
            if quiet_dd.value is None or quiet_dd.value == "":
                print("No quiet FITS selected.")
                return

            quiet_f = str(quiet_dd.value)
            roi = dict(
                xmin=float(xmin_w.value),
                xmax=float(xmax_w.value),
                ymin=float(ymin_w.value),
                ymax=float(ymax_w.value),
            )

            crop_hw = float(crop_in.value)
            clim = (float(clim_low.value), float(clim_high.value))

            print("Quiet FITS:", os.path.basename(quiet_f))
            print("ROI (arcsec):", roi)
            print("Crop half-width (as):", crop_hw)
            print("CLim pct:", clim)
            print(f"Fit params: thresh_frac={thresh_frac}, min_points={min_points}")

            # preview + centroid
            m_sub, fit, obstime, freq = _plot_quiet_preview(quiet_f, roi, crop_hw, clim)

            if fit is None:
                print("\n[FAIL] centroid fit failed. Try adjust ROI or thresh/min_points.")
                return

            # dx/dy
            dx_pix, dy_pix, cen_pix, tgt_pix = _solve_shift(m_sub, fit)
            print("\ncentroid (arcsec):", fit["cen_Tx_as"], fit["cen_Ty_as"])
            print("centroid pix(full cropped):", cen_pix)
            print("target pix(full cropped):", tgt_pix)
            print("dx_pix, dy_pix:", dx_pix, dy_pix)

            # save json
            tag = _sb_tag(sb_dd.value)
            sb_out = _ensure_dir(out_root / tag)
            sol_path = sb_out / "step5a_solution.json"

            sol = dict(
                step="5A",
                created_utc=datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
                step4_root=str(step4_root),
                sb=tag,
                quiet_fits=os.path.basename(quiet_f),
                quiet_fits_path=str(Path(quiet_f).resolve()),
                crop_half_width_arcsec=crop_hw,
                roi_arcsec=roi,
                thresh_frac=float(thresh_frac),
                min_points=int(min_points),
                centroid_arcsec=dict(Tx=fit["cen_Tx_as"], Ty=fit["cen_Ty_as"]),
                dx_pix=float(dx_pix),
                dy_pix=float(dy_pix),
                note="dx_pix,dy_pix defined on CROPPED+ROTATED HPC map grid; apply as scipy.ndimage.shift(data, shift=(dy,dx))",
            )

            with open(sol_path, "w") as f:
                json.dump(sol, f, indent=2)

            print("\nSaved:", sol_path)

    refresh_btn.on_click(_refresh)
    sb_dd.observe(lambda ch: _refresh(), names="value")
    run_btn.on_click(_on_run)

    _refresh()

    display(
        w.VBox([
            w.HBox([sb_dd, refresh_btn]),
            quiet_dd,
            w.HBox([crop_in, clim_low, clim_high]),
            w.HBox([xmin_w, xmax_w, ymin_w, ymax_w]),
            run_btn,
            out
        ])
    )

def run_step5b_iocorrect_apply_ui(
    step4_root,
    out_root=None,
    default_sb="SB410",
    default_crop_half_width_arcsec=5000,
    default_clim_pct=(5, 99),
    default_contours=(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95),
):
    """
    Step-5B (IO-correct APPLY, WCS mode):
      - read step5a_solution.json (per SB) -> centroid offset (Tx,Ty) and dx/dy
      - apply correction by adjusting WSClean FITS WCS (CRVAL1/2) ONLY (no data shifting)
      - write corrected FITS to out_root/SBxxx/corr_fits/
      - generate step4-style quicklook BUT side-by-side: BEFORE vs AFTER(WCS corrected)
      - optional video mp4/gif
    """

    import os, re, glob, json, shutil, subprocess
    import numpy as np
    import ipywidgets as w
    from pathlib import Path
    from IPython.display import display, clear_output

    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    import astropy.units as u
    from astropy.io import fits
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation

    import sunpy.map
    from sunpy.coordinates import frames, sun

    # ---------------- helpers ----------------
    def _sb_tag(sb):
        return sb.replace(".MS", "").replace("SB", "SB")

    def _find_sbs(step4_root_path: Path):
        sbs = []
        for p in sorted(step4_root_path.glob("SB*")):
            if p.is_dir() and re.match(r"SB\d{3}$", p.name):
                sbs.append(p.name + ".MS")
        return sbs

    def _list_fits(step4_root_path: Path, sb: str):
        tag = _sb_tag(sb)  # "SB410"
        sb_dir = step4_root_path / tag
        pats = [
            str(sb_dir / f"{tag}-t*-image.fits"),
            str(sb_dir / f"{tag}_t*-image.fits"),
            str(sb_dir / "*-image.fits"),
        ]
        files = []
        for pat in pats:
            files.extend(glob.glob(pat))
        files = sorted(set(files))
        return files

    def _parse_tindex(fname: str):
        m = re.search(r"-t(\d+)-", os.path.basename(fname))
        return int(m.group(1)) if m else 10**9

    def _load_sort_fits(files):
        return sorted(files, key=lambda x: (_parse_tindex(x), os.path.basename(x)))

    def _ensure_dir(p: Path):
        p.mkdir(parents=True, exist_ok=True)
        return p

    def _parse_contours(s: str):
        s = (s or "").strip()
        if not s:
            return []
        vals = []
        for x in s.split(","):
            x = x.strip()
            if not x:
                continue
            try:
                vals.append(float(x))
            except Exception:
                pass
        return vals

    def _guess_obstime(hdr):
        for k in ["DATE-OBS", "DATEOBS", "DATE_OBS", "DATE"]:
            if k in hdr:
                try:
                    return Time(hdr[k])
                except Exception:
                    pass
        return Time.now()

    # ----- core: build HPC rotated cropped map exactly like step4 quicklook -----
    def _build_hpc_rot_submap_from_datahdr(data2d, hdr, crop_half_width_arcsec, site):
        # ensure 2D
        data = data2d
        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image, shape={getattr(data,'shape',None)}")

        obstime = _guess_obstime(hdr)

        freq_Hz = hdr.get("CRVAL3", None)
        frequency = (freq_Hz * u.Hz) if freq_Hz is not None else (np.nan * u.Hz)

        # pixel scale
        cdelt1 = abs(hdr.get("CDELT1", np.nan)) * u.deg
        cdelt2 = abs(hdr.get("CDELT2", np.nan)) * u.deg
        cdelt1 = cdelt1.to(u.arcsec) if np.isfinite(cdelt1.value) else (np.nan * u.arcsec)
        cdelt2 = cdelt2.to(u.arcsec) if np.isfinite(cdelt2.value) else (np.nan * u.arcsec)

        # observer
        site_gcrs = SkyCoord(site.get_gcrs(obstime))

        cunit1 = u.Unit(hdr.get("CUNIT1", "deg"))
        cunit2 = u.Unit(hdr.get("CUNIT2", "deg"))

        # reference direction from CRVAL1/2
        ref_gcrs = SkyCoord(
            hdr["CRVAL1"] * cunit1,
            hdr["CRVAL2"] * cunit2,
            frame="gcrs",
            obstime=obstime,
            obsgeoloc=site_gcrs.cartesian,
            obsgeovel=site_gcrs.velocity.to_cartesian(),
            distance=site_gcrs.hcrs.distance,
        )

        ref_hpc = ref_gcrs.transform_to(frames.Helioprojective(observer=site_gcrs))

        # rotate north-up
        P1 = sun.P(obstime)

        ref_pix = np.array([hdr["CRPIX1"] - 1, hdr["CRPIX2"] - 1]) * u.pixel
        scale = np.array([cdelt1.value, cdelt2.value]) * (u.arcsec / u.pixel)

        new_header = sunpy.map.make_fitswcs_header(
            data=data,
            coordinate=ref_hpc,
            reference_pixel=ref_pix,
            scale=scale,
            rotation_angle=-P1,
            wavelength=frequency.to(u.MHz) if np.isfinite(frequency.value) else None,
            observatory="NenuFAR (Nançay)",
        )

        rmap = sunpy.map.Map(data, new_header)
        rmap_rot = rmap.rotate()

        # crop
        hw = float(crop_half_width_arcsec) * u.arcsec
        bl = SkyCoord(-hw, -hw, frame=rmap_rot.coordinate_frame)
        tr = SkyCoord(hw, hw, frame=rmap_rot.coordinate_frame)
        try:
            rsub = rmap_rot.submap(bl, top_right=tr)
        except Exception:
            rsub = rmap_rot

        return rsub, obstime, frequency

    # ----- WCS correction (B): update CRVAL1/2 by shifting ref_hpc by (-Tx,-Ty) then back to GCRS -----
    def _apply_wcs_correction_to_header(hdr_in, Tx_as, Ty_as, site):
        hdr = hdr_in.copy()
        obstime = _guess_obstime(hdr)

        site_gcrs = SkyCoord(site.get_gcrs(obstime))

        cunit1 = u.Unit(hdr.get("CUNIT1", "deg"))
        cunit2 = u.Unit(hdr.get("CUNIT2", "deg"))

        ref_gcrs_old = SkyCoord(
            hdr["CRVAL1"] * cunit1,
            hdr["CRVAL2"] * cunit2,
            frame="gcrs",
            obstime=obstime,
            obsgeoloc=site_gcrs.cartesian,
            obsgeovel=site_gcrs.velocity.to_cartesian(),
            distance=site_gcrs.hcrs.distance,
        )
        ref_hpc_old = ref_gcrs_old.transform_to(frames.Helioprojective(observer=site_gcrs))

        # shift reference in HPC tangent plane
        ref_hpc_new = SkyCoord(
            (ref_hpc_old.Tx - Tx_as * u.arcsec),
            (ref_hpc_old.Ty - Ty_as * u.arcsec),
            frame=ref_hpc_old.frame,
        )
        ref_gcrs_new = ref_hpc_new.transform_to(ref_gcrs_old.frame)

        # write back CRVAL1/2 in original units
        hdr["CRVAL1"] = ref_gcrs_new.ra.to_value(cunit1)
        hdr["CRVAL2"] = ref_gcrs_new.dec.to_value(cunit2)

        hdr.add_history(f"STEP5B WCS-correct: CRVAL shifted via HPC offset Tx={Tx_as:.6f} arcsec, Ty={Ty_as:.6f} arcsec")
        return hdr

    # ----- plotting: BEFORE vs AFTER (both step4 style) -----
    def _quicklook_before_after(
        fits_path,
        out_png,
        crop_half_width_arcsec,
        clim_pct,
        contours_perc,
        site,
        overwrite=False,
        cmap="viridis",
        draw_beam=True,
        hdr_after=None,   # if provided, use (data, hdr_after) for AFTER
    ):
        out_png = str(out_png)
        if (not overwrite) and os.path.exists(out_png):
            return out_png

        with fits.open(fits_path) as hdul:
            hdr0 = hdul[0].header
            data0 = hdul[0].data

        # ensure 2D for both
        data0_2d = data0
        while data0_2d is not None and getattr(data0_2d, "ndim", 0) > 2:
            data0_2d = data0_2d[0]
        data0_2d = np.squeeze(data0_2d)

        # BEFORE map
        m_before, obstime, frequency = _build_hpc_rot_submap_from_datahdr(
            data0_2d, hdr0, crop_half_width_arcsec, site
        )

        # AFTER map (same data, corrected header)
        hdr1 = hdr_after if hdr_after is not None else hdr0
        m_after, _, _ = _build_hpc_rot_submap_from_datahdr(
            data0_2d, hdr1, crop_half_width_arcsec, site
        )

        # shared scaling from BEFORE
        vmin, vmax = np.nanpercentile(m_before.data, clim_pct)

        fig = plt.figure(figsize=(14.5, 7.0), constrained_layout=False)
        gs = fig.add_gridspec(nrows=1, ncols=3, width_ratios=[1.0, 1.0, 0.04], wspace=0.15)

        axL = fig.add_subplot(gs[0, 0], projection=m_before)
        axR = fig.add_subplot(gs[0, 1], projection=m_after)
        cax = fig.add_subplot(gs[0, 2])

        imL = m_before.plot(axes=axL, cmap=cmap, vmin=vmin, vmax=vmax)
        imR = m_after.plot(axes=axR, cmap=cmap, vmin=vmin, vmax=vmax)

        for ax, mm in [(axL, m_before), (axR, m_after)]:
            try:
                mm.draw_limb(axes=ax)
            except Exception:
                pass
            try:
                mm.draw_grid(axes=ax)
            except Exception:
                pass
            if contours_perc:
                try:
                    mm.draw_contours(np.array(contours_perc) * u.percent, colors="k",
                                     linewidths=1.1, alpha=0.8, axes=ax)
                except Exception:
                    pass

        axL.set_title("BEFORE (step4)", fontsize=13)
        axR.set_title("AFTER (WCS-corrected)", fontsize=13)

        # colorbar
        cbar = fig.colorbar(imL, cax=cax)
        sf = ScalarFormatter(useMathText=True)
        sf.set_powerlimits((-2, 3))
        cbar.formatter = sf
        cbar.update_ticks()
        cbar.set_label("Stokes I (arb.)", labelpad=12)

        # titles
        left_title = f"{os.path.basename(fits_path)}   {obstime.isot}"
        fig.text(0.02, 0.98, left_title, ha="left", va="top", fontsize=16)
        if np.isfinite(frequency.value):
            freq_str = f"{frequency.to_value(u.MHz):.1f} MHz"
            fig.text(0.98, 0.98, freq_str, ha="right", va="top",
                     fontsize=15, bbox=dict(fc=(1, 1, 1, 0.7), ec="0.7", pad=2))

        # optional beam on LEFT only (same as step4)
        if draw_beam:
            try:
                from matplotlib.patches import Ellipse
                bmaj = hdr0.get("BMAJ"); bmin = hdr0.get("BMIN"); bpa = hdr0.get("BPA")
                if bmaj is not None and bmin is not None and bpa is not None:
                    bmaj_as = (abs(bmaj) * u.deg).to_value(u.arcsec)
                    bmin_as = (abs(bmin) * u.deg).to_value(u.arcsec)
                    x0 = axL.get_xlim()[0] + 0.08*(axL.get_xlim()[1]-axL.get_xlim()[0])
                    y0 = axL.get_ylim()[0] + 0.08*(axL.get_ylim()[1]-axL.get_ylim()[0])
                    e = Ellipse((x0, y0), width=bmaj_as, height=bmin_as, angle=-float(bpa),
                                edgecolor="w", facecolor="none", lw=1.5)
                    axL.add_patch(e)
            except Exception:
                pass

        Path(out_png).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_png

    def _make_video_ffmpeg(png_list, out_mp4, fps=8):
        png_list = list(png_list)
        if len(png_list) == 0:
            return None
        out_mp4 = Path(out_mp4)

        frames_dir = out_mp4.parent / (out_mp4.stem + "_frames")
        if frames_dir.exists():
            shutil.rmtree(frames_dir)
        frames_dir.mkdir(parents=True, exist_ok=True)

        png_list = sorted(png_list)
        for i, p in enumerate(png_list):
            p = Path(p)
            link_path = frames_dir / f"{i:06d}.png"
            try:
                os.symlink(p, link_path)
            except Exception:
                shutil.copy2(p, link_path)

        cmd = [
            "ffmpeg", "-y",
            "-framerate", str(int(fps)),
            "-i", str(frames_dir / "%06d.png"),
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-movflags", "+faststart",
            str(out_mp4),
        ]
        subprocess.check_call(cmd)
        return str(out_mp4)

    def _make_video_imageio(png_list, out_gif, fps=8):
        import imageio.v2 as imageio
        imgs = [imageio.imread(p) for p in png_list]
        out_gif = str(out_gif)
        imageio.mimsave(out_gif, imgs, duration=1.0/max(fps, 1))
        return out_gif

    # ---------------- paths ----------------
# ---------------- paths ----------------
    step4_root = Path(step4_root)

    if out_root is None:
        # 从 step4_root 里提取事件日期，例如 20240310
        import re
        m = re.search(r"(20\d{6})", str(step4_root))
        if m:
            event_tag = m.group(1)
        else:
            # fallback（极端情况）
            import datetime
            event_tag = datetime.datetime.now().strftime("%Y%m%d")

        out_root = Path("/data/jzhang/nenufar_workflows") / f"step_iocorrect_outputs_{event_tag}"

    out_root = Path(out_root)

    sbs = _find_sbs(step4_root)
    if len(sbs) == 0:
        raise ValueError(f"No SBxxx directories found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    # site fixed to NenuFAR (same as step4)
    site = EarthLocation(lat=47.382*u.deg, lon=2.195*u.deg, height=136.0*u.m)

    # ---------------- UI widgets ----------------
    sb_dd = w.Dropdown(options=sbs, value=default_sb_ms, description="SB")
    refresh_btn = w.Button(description="Refresh FITS", button_style="")
    fits_sel = w.SelectMultiple(options=[], description="FITS", layout=w.Layout(width="95%", height="220px"))

    select_all_cb = w.Checkbox(value=True, description="Select ALL image.fits (ignore manual selection)")
    overwrite_cb = w.Checkbox(value=False, description="Overwrite outputs")
    make_video_cb = w.Checkbox(value=False, description="Make video (mp4/gif)")
    fps_in = w.IntText(value=8, description="FPS")

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop hw (as)")
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low%")
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high%")
    contours_in = w.Text(value=",".join(map(str, default_contours)), description="Contours%")

    run_btn = w.Button(description="Run Step-5B (apply WCS)", button_style="success")
    out = w.Output()

    # ---------------- behaviors ----------------
    def _refresh(_=None):
        files = _list_fits(step4_root, sb_dd.value)
        files = _load_sort_fits(files)
        fits_sel.options = files
        with out:
            out.clear_output()
            print("Step4 root:", step4_root)
            print("SB:", sb_dd.value, "->", len(files), "image.fits")
            if len(files) > 0:
                print("First:", os.path.basename(files[0]))
                print("Last :", os.path.basename(files[-1]))

    def _load_solution_json(sb_tag: str):
        # expected location from your step5a:
        # out_root/SBxxx/step5a_solution.json
        sol = out_root / sb_tag / "step5a_solution.json"
        if not sol.exists():
            raise FileNotFoundError(f"step5a_solution.json not found: {sol}")
        with open(sol, "r") as f:
            return json.load(f), sol

    def _on_run(_):
        with out:
            out.clear_output()

            sb = sb_dd.value
            tag = _sb_tag(sb)

            # fits selection
            all_files = list(fits_sel.options)
            chosen = list(all_files) if select_all_cb.value else list(fits_sel.value)
            if len(chosen) == 0:
                print("No FITS selected. Either tick 'Select ALL' or select manually.")
                return

            # load step5a solution
            try:
                sol, sol_path = _load_solution_json(tag)
            except Exception as e:
                print("[FAIL] cannot load step5a solution:", e)
                return

            # ---------- robust centroid parser ----------
            def _as_float(x):
                try:
                    if x is None:
                        return np.nan
                    if isinstance(x, (np.floating, float, int, np.integer)):
                        return float(x)
                    return float(str(x).strip())
                except Exception:
                    return np.nan

            def _parse_centroid(sol_dict):
                # Priority 1: explicit keys
                Tx = _as_float(sol_dict.get("cen_Tx_as", np.nan))
                Ty = _as_float(sol_dict.get("cen_Ty_as", np.nan))
                if np.isfinite(Tx) and np.isfinite(Ty):
                    return Tx, Ty

                # Priority 2: centroid_arcsec field (can be list/tuple/dict/str)
                if "centroid_arcsec" in sol_dict:
                    c = sol_dict["centroid_arcsec"]

                    # list/tuple/ndarray
                    if isinstance(c, (list, tuple, np.ndarray)) and len(c) >= 2:
                        Tx = _as_float(c[0])
                        Ty = _as_float(c[1])
                        if np.isfinite(Tx) and np.isfinite(Ty):
                            return Tx, Ty

                    # dict with common keys
                    if isinstance(c, dict):
                        for kx, ky in [("Tx", "Ty"), ("tx", "ty"), ("x", "y"), ("cen_Tx_as", "cen_Ty_as")]:
                            Tx = _as_float(c.get(kx, np.nan))
                            Ty = _as_float(c.get(ky, np.nan))
                            if np.isfinite(Tx) and np.isfinite(Ty):
                                return Tx, Ty

                    # string like "(-441.9, 699.0)" or "-441.9,699.0"
                    if isinstance(c, str):
                        s = c.strip().replace("(", "").replace(")", "")
                        parts = [p.strip() for p in s.split(",")]
                        if len(parts) >= 2:
                            Tx = _as_float(parts[0])
                            Ty = _as_float(parts[1])
                            if np.isfinite(Tx) and np.isfinite(Ty):
                                return Tx, Ty

                # Priority 3: other possible names (future-proof)
                for kx, ky in [("Tx_as", "Ty_as"), ("Tx", "Ty"), ("centroid_Tx_as", "centroid_Ty_as")]:
                    Tx = _as_float(sol_dict.get(kx, np.nan))
                    Ty = _as_float(sol_dict.get(ky, np.nan))
                    if np.isfinite(Tx) and np.isfinite(Ty):
                        return Tx, Ty

                return np.nan, np.nan

            Tx_as, Ty_as = _parse_centroid(sol)
            if (not np.isfinite(Tx_as)) or (not np.isfinite(Ty_as)):
                print("[FAIL] cannot parse centroid from solution JSON:", sol_path)
                print("Available keys:", list(sol.keys()))
                return

            # Your saved format from step5a screenshot:
            # - centroid (arcsec) printed
            # - dx_pix/dy_pix printed
            # But JSON may store keys slightly differently; we tolerate both.
            if not np.isfinite(Tx_as) or not np.isfinite(Ty_as):
                # try older keys
                Tx_as = float(sol.get("cen_Tx_as", np.nan))
                Ty_as = float(sol.get("cen_Ty_as", np.nan))

            if not np.isfinite(Tx_as) or not np.isfinite(Ty_as):
                print("[FAIL] centroid arcsec not found in solution JSON:", sol_path)
                print("keys:", list(sol.keys()))
                return

            # params
            clim = (float(clim_low.value), float(clim_high.value))
            contours = _parse_contours(contours_in.value)
            crop_hw = float(crop_in.value)

            # outputs
            sb_out = _ensure_dir(out_root / tag)
            corr_dir = _ensure_dir(sb_out / "corr_fits")
            ql_dir = _ensure_dir(sb_out / "quicklook_step5b")
            log_path = sb_out / "step5b_apply.log"

            print("SB:", sb)
            print("Using solution:", sol_path)
            print(f"centroid offset (arcsec): Tx={Tx_as:.3f}, Ty={Ty_as:.3f}")
            print("Selected FITS:", len(chosen), "(select_all =", select_all_cb.value, ")")
            print("Out:", sb_out)
            print("corr_fits:", corr_dir)
            print("quicklooks:", ql_dir)
            print("Overwrite:", overwrite_cb.value, "| video:", make_video_cb.value, "fps:", fps_in.value)

            pngs = []
            ok = 0

            for fpath in chosen:
                fpath = str(fpath)
                fin = Path(fpath)
                fout = corr_dir / fin.name.replace("-image.fits", "-image_corrWCS.fits")

                # write corrected FITS (WCS only)
                try:
                    if (not overwrite_cb.value) and fout.exists():
                        pass
                    else:
                        with fits.open(fin) as hdul:
                            hdr0 = hdul[0].header
                            data0 = hdul[0].data
                        hdr1 = _apply_wcs_correction_to_header(hdr0, Tx_as=Tx_as, Ty_as=Ty_as, site=site)

                        hdu = fits.PrimaryHDU(data=data0, header=hdr1)
                        hdul_out = fits.HDUList([hdu])
                        hdul_out.writeto(fout, overwrite=True)

                    # quicklook before/after side-by-side
                    out_png = ql_dir / (fin.name.replace("-image.fits", "_step5b_before_after.png"))
                    with fits.open(fin) as hdul:
                        hdr0 = hdul[0].header
                        data0 = hdul[0].data
                    hdr_after = _apply_wcs_correction_to_header(hdr0, Tx_as=Tx_as, Ty_as=Ty_as, site=site)

                    png = _quicklook_before_after(
                        fits_path=str(fin),
                        out_png=str(out_png),
                        crop_half_width_arcsec=crop_hw,
                        clim_pct=clim,
                        contours_perc=contours,
                        site=site,
                        overwrite=overwrite_cb.value,
                        hdr_after=hdr_after,
                    )
                    pngs.append(png)
                    ok += 1
                    with open(log_path, "a") as f:
                        f.write(f"OK  {fin} -> {fout} | {png}\n")
                except Exception as e:
                    with open(log_path, "a") as f:
                        f.write(f"FAIL {fin} : {e}\n")
                    print("[FAIL]", fin.name, ":", e)

            print(f"Done. OK {ok}/{len(chosen)}")
            print("Log:", log_path)

            # video
            if make_video_cb.value and len(pngs) > 1:
                mp4 = sb_out / f"{tag}_step5b_before_after.mp4"
                gif = sb_out / f"{tag}_step5b_before_after.gif"
                try:
                    if shutil.which("ffmpeg"):
                        _make_video_ffmpeg(pngs, mp4, fps=int(fps_in.value))
                        print("Video:", mp4)
                    else:
                        raise RuntimeError("ffmpeg not found")
                except Exception:
                    try:
                        _make_video_imageio(pngs, gif, fps=int(fps_in.value))
                        print("GIF:", gif)
                    except Exception as e:
                        print("[WARN] video failed:", e)

    refresh_btn.on_click(_refresh)
    run_btn.on_click(_on_run)

    _refresh()

    display(
        w.VBox([
            w.HBox([sb_dd, refresh_btn, select_all_cb, overwrite_cb, make_video_cb, fps_in]),
            fits_sel,
            w.HBox([crop_in, clim_low, clim_high]),
            w.HBox([contours_in]),
            run_btn,
            out
        ])
    )

def run_step5c_centroid_ui(
    step4_root,
    step5b_root=None,   # optional: /data/.../step_iocorrect_outputs_20240310
    out_root=None,      # default: /data/.../step_iocentroid_outputs_YYYYMMDD
    default_sb="SB410",
    default_crop_half_width_arcsec=5000,
    default_clim_pct=(5, 99),
    default_contours=(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
    default_roi=(-1500, 1000, 0, 1500),  # xmin,xmax,ymin,ymax arcsec
    default_thresh_frac=0.5,
    default_min_points=30,
    site_lat_deg=47.382,
    site_lon_deg=2.195,
    site_height_m=136.0,
):
    """
    Step-5C Centroid UI:
      - choose FITS (raw step4 image.fits OR step5b corrected corr_fits/*.fits)
      - define ROI and fit centroid via 2D Gaussian on pixels above thresh_frac*peak
      - generate step4-style quicklooks and save results table
      - optional: make movie (mp4/gif) from generated pngs
    """
    import os
    import re
    import json
    import glob
    import shutil
    import subprocess
    import datetime
    from pathlib import Path

    import numpy as np
    import ipywidgets as w
    from IPython.display import display, clear_output

    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.io import fits
    from astropy.time import Time

    import sunpy.map
    from sunpy.coordinates import frames, sun

    # scipy only for gaussian fit
    from scipy.optimize import curve_fit

    # -------------------------
    # helpers: consistent paths
    # -------------------------
    def _infer_event_tag_from_step4_root(step4_root_path: Path) -> str:
        s = str(step4_root_path)
        m = re.search(r"(20\d{6})", s)
        if m:
            return m.group(1)
        return datetime.datetime.now().strftime("%Y%m%d")

    def _default_out_root(prefix: str, event_tag: str) -> Path:
        return Path("/data/jzhang/nenufar_workflows") / f"{prefix}_{event_tag}"

    # -------------------------
    # helpers: SB + FITS listing
    # -------------------------
    def _sb_tag(sb):
        return sb.replace(".MS", "").replace("SB", "SB")

    def _find_sbs(root: Path):
        sbs = []
        for p in sorted(root.glob("SB*")):
            if p.is_dir() and re.match(r"SB\d{3}$", p.name):
                sbs.append(p.name + ".MS")
        return sbs

    def _parse_tindex(fname: str):
        m = re.search(r"-t(\d+)-", os.path.basename(fname))
        return int(m.group(1)) if m else 10**9

    def _load_sort_fits(files):
        return sorted(files, key=lambda x: (_parse_tindex(x), os.path.basename(x)))

    def _list_step4_image_fits(step4_root_path: Path, sb: str):
        """ONLY list *image.fits (no dirty/psf/model/residual)."""
        tag = _sb_tag(sb)
        sb_dir = step4_root_path / tag
        pats = [
            str(sb_dir / f"{tag}-t*-image.fits"),
            str(sb_dir / f"{tag}_t*-image.fits"),
            str(sb_dir / "*-image.fits"),
        ]
        files = []
        for pat in pats:
            files.extend(glob.glob(pat))
        files = sorted(set(files))
        return [f for f in files if os.path.isfile(f)]

    def _list_step5b_corr_fits(step5b_root_path: Path, sb: str):
        if step5b_root_path is None:
            return []
        tag = _sb_tag(sb)
        corr_dir = step5b_root_path / tag / "corr_fits"
        if not corr_dir.exists():
            return []
        files = sorted([str(p) for p in corr_dir.glob("*.fits") if p.is_file()])
        return files

    def _ensure_dir(p: Path):
        p.mkdir(parents=True, exist_ok=True)
        return p

    # -------------------------
    # core: FITS -> HPC SunPy map (rotated north-up) + crop
    # -------------------------
    site = EarthLocation(lat=float(site_lat_deg)*u.deg, lon=float(site_lon_deg)*u.deg, height=float(site_height_m)*u.m)

    def _read_2d(fits_path: str):
        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image: {fits_path} shape={getattr(data,'shape',None)}")
        return data.astype(float), hdr

    def _get_obstime(hdr):
        for k in ["DATE-OBS", "DATEOBS", "DATE_OBS", "DATE"]:
            if k in hdr:
                try:
                    return Time(hdr[k])
                except Exception:
                    pass
        return Time.now()

    def _build_hpc_rotated_map(data, hdr):
        obstime = _get_obstime(hdr)
        freq_Hz = hdr.get("CRVAL3", np.nan)
        frequency = (freq_Hz * u.Hz) if np.isfinite(freq_Hz) else (np.nan * u.Hz)

        site_gcrs = SkyCoord(site.get_gcrs(obstime))

        cunit1 = u.Unit(hdr.get("CUNIT1", "deg"))
        cunit2 = u.Unit(hdr.get("CUNIT2", "deg"))

        ref_gcrs = SkyCoord(
            hdr["CRVAL1"] * cunit1,
            hdr["CRVAL2"] * cunit2,
            frame="gcrs",
            obstime=obstime,
            obsgeoloc=site_gcrs.cartesian,
            obsgeovel=site_gcrs.velocity.to_cartesian(),
            distance=site_gcrs.hcrs.distance,
        )

        ref_hpc = ref_gcrs.transform_to(frames.Helioprojective(observer=site_gcrs))

        cdelt1 = (abs(hdr.get("CDELT1", np.nan)) * u.deg).to(u.arcsec)
        cdelt2 = (abs(hdr.get("CDELT2", np.nan)) * u.deg).to(u.arcsec)

        P1 = sun.P(obstime)

        ref_pix = np.array([hdr["CRPIX1"] - 1, hdr["CRPIX2"] - 1]) * u.pixel
        scale = np.array([cdelt1.value, cdelt2.value]) * (u.arcsec / u.pixel)

        new_header = sunpy.map.make_fitswcs_header(
            data=data,
            coordinate=ref_hpc,
            reference_pixel=ref_pix,
            scale=scale,
            rotation_angle=-P1,
            wavelength=frequency.to(u.MHz) if np.isfinite(frequency.value) else None,
            observatory="NenuFAR (Nançay)",
        )

        m = sunpy.map.Map(data, new_header)
        m = m.rotate()
        return m, obstime, frequency.to(u.MHz) if np.isfinite(frequency.value) else (np.nan * u.MHz)

    def _crop_center(m, crop_hw_as):
        hw = float(crop_hw_as) * u.arcsec
        bl = SkyCoord(-hw, -hw, frame=m.coordinate_frame)
        tr = SkyCoord(+hw, +hw, frame=m.coordinate_frame)
        try:
            return m.submap(bl, top_right=tr)
        except Exception:
            return m

    def _submap_from_roi(m, roi):
        xmin, xmax, ymin, ymax = roi
        p0 = SkyCoord(float(xmin)*u.arcsec, float(ymin)*u.arcsec, frame=m.coordinate_frame)
        p1 = SkyCoord(float(xmax)*u.arcsec, float(ymax)*u.arcsec, frame=m.coordinate_frame)
        return m.submap(p0, top_right=p1)

    # -------------------------
    # centroid fit: 2D Gaussian (rotated)
    # -------------------------
    def _centroid_gauss(roi_map, thresh_frac=0.5, min_points=30):
        z2 = np.array(roi_map.data, dtype=float)
        if np.isnan(z2).any():
            finite = z2[np.isfinite(z2)]
            z2[np.isnan(z2)] = np.nanmin(finite) if finite.size else 0.0

        amp0 = np.nanmax(z2)
        if not np.isfinite(amp0) or amp0 <= 0:
            return None

        thr = float(thresh_frac) * amp0
        mask = z2 > thr
        if np.count_nonzero(mask) < int(min_points):
            return None

        yy, xx = np.indices(z2.shape)
        x = xx[mask].ravel().astype(float)
        y = yy[mask].ravel().astype(float)
        z = z2[mask].ravel().astype(float)

        def gauss2d(coords, A, x0, y0, theta, sx, sy):
            x_, y_ = coords
            ct, st = np.cos(theta), np.sin(theta)
            xp = ct*(x_ - x0) + st*(y_ - y0)
            yp = -st*(x_ - x0) + ct*(y_ - y0)
            return (A * np.exp(-0.5*((xp/sx)**2 + (yp/sy)**2))).ravel()

        y0_idx, x0_idx = np.unravel_index(np.nanargmax(z2), z2.shape)
        x0g, y0g = float(x0_idx), float(y0_idx)

        s_guess = max(2.0, 0.15 * min(z2.shape))
        p0 = [amp0, x0g, y0g, 0.0, s_guess, s_guess]
        bounds = (
            [0.0, 0.0, 0.0, -np.pi, 1.0, 1.0],
            [3*amp0, z2.shape[1]-1, z2.shape[0]-1, np.pi, max(z2.shape), max(z2.shape)],
        )

        try:
            popt, pcov = curve_fit(gauss2d, (x, y), z, p0=p0, bounds=bounds, maxfev=20000)
        except Exception:
            return None

        A, x0_fit, y0_fit, theta_fit, sx_fit, sy_fit = popt
        cen_world = roi_map.pixel_to_world(x0_fit * u.pix, y0_fit * u.pix)

        # Tx/Ty arcsec
        try:
            Tx_as = float(cen_world.Tx.to_value(u.arcsec))
            Ty_as = float(cen_world.Ty.to_value(u.arcsec))
        except Exception:
            Tx_as = float(cen_world.spherical.lon.to_value(u.arcsec))
            Ty_as = float(cen_world.spherical.lat.to_value(u.arcsec))

        return dict(
            cen_world=cen_world,
            cen_Tx_as=Tx_as,
            cen_Ty_as=Ty_as,
            popt=[float(v) for v in popt],
        )

    # -------------------------
    # plotting: step4-like + ROI box + contours + centroid marker
    # -------------------------
    def _parse_contours(s: str):
        s = (s or "").strip()
        if not s:
            return []
        vals = []
        for x in s.split(","):
            x = x.strip()
            if not x:
                continue
            try:
                vals.append(float(x))
            except Exception:
                pass
        return vals

    def _draw_roi_box(ax, m, roi):
        xmin, xmax, ymin, ymax = roi
        p0 = SkyCoord(float(xmin)*u.arcsec, float(ymin)*u.arcsec, frame=m.coordinate_frame)
        p1 = SkyCoord(float(xmax)*u.arcsec, float(ymax)*u.arcsec, frame=m.coordinate_frame)
        x0, y0 = m.world_to_pixel(p0)
        x1, y1 = m.world_to_pixel(p1)
        x0 = float(np.atleast_1d(getattr(x0, "value", x0))[0])
        y0 = float(np.atleast_1d(getattr(y0, "value", y0))[0])
        x1 = float(np.atleast_1d(getattr(x1, "value", x1))[0])
        y1 = float(np.atleast_1d(getattr(y1, "value", y1))[0])
        ax.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0], transform=ax.get_transform("pixel"), lw=2)

    def _quicklook_centroid_one(fits_path, out_png, crop_hw_as, clim_pct, contours, roi, thresh_frac, min_points, overwrite=False):
        if (not overwrite) and os.path.exists(out_png):
            return out_png, None

        data, hdr = _read_2d(fits_path)
        m0, obstime, freq_mhz = _build_hpc_rotated_map(data, hdr)
        m = _crop_center(m0, crop_hw_as)
        roi_m = _submap_from_roi(m, roi)

        fit = _centroid_gauss(roi_m, thresh_frac=thresh_frac, min_points=min_points)

        vmin, vmax = np.nanpercentile(m.data, clim_pct)

        fig = plt.figure(figsize=(8.8, 7.4), constrained_layout=False)
        gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[1.0, 0.04], wspace=0.15)
        ax = fig.add_subplot(gs[0, 0], projection=m)
        cax = fig.add_subplot(gs[0, 1])

        im = m.plot(axes=ax, cmap="viridis", vmin=vmin, vmax=vmax)
        try:
            m.draw_limb(axes=ax)
            m.draw_grid(axes=ax)
        except Exception:
            pass

        if contours:
            try:
                m.draw_contours(np.array(contours) * u.percent, colors="k", linewidths=1.1, alpha=0.8, axes=ax)
            except Exception:
                pass

        _draw_roi_box(ax, m, roi)

        # disk center marker
        c0 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=m.coordinate_frame)
        ax.plot_coord(c0, marker="x", markersize=10, mew=2)

        # centroid marker
        if fit is not None:
            ax.plot_coord(fit["cen_world"], marker="+", markersize=14, mew=2)
            ax.text(
                0.02, 0.02,
                f'cen=({fit["cen_Tx_as"]:.1f}", {fit["cen_Ty_as"]:.1f}")',
                transform=ax.transAxes, ha="left", va="bottom",
                bbox=dict(fc=(1,1,1,0.6), ec="none", pad=2),
            )

        # colorbar
        cbar = fig.colorbar(im, cax=cax)
        sf = ScalarFormatter(useMathText=True)
        sf.set_powerlimits((-2, 3))
        cbar.formatter = sf
        cbar.update_ticks()
        cbar.set_label("Stokes I (arb.)", labelpad=12)

        # titles
        left_title = f"{os.path.basename(fits_path)}   {obstime.isot}"
        fig.text(0.02, 0.98, left_title, ha="left", va="top", fontsize=16)
        if np.isfinite(freq_mhz.value):
            fig.text(0.98, 0.98, f"{freq_mhz.to_value(u.MHz):.1f} MHz",
                     ha="right", va="top", fontsize=15,
                     bbox=dict(fc=(1, 1, 1, 0.7), ec="0.7", pad=2))

        _ensure_dir(Path(out_png).parent)
        plt.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)

        return out_png, dict(
            obstime_isot=str(obstime.isot),
            freq_mhz=float(freq_mhz.to_value(u.MHz)) if np.isfinite(freq_mhz.value) else np.nan,
            cen_Tx_as=(fit["cen_Tx_as"] if fit else np.nan),
            cen_Ty_as=(fit["cen_Ty_as"] if fit else np.nan),
            ok=(fit is not None),
        )

    # -------------------------
    # movie makers (same spirit as step4/5b)
    # -------------------------
    def _make_video_ffmpeg(png_list, out_mp4, fps=8):
        png_list = list(png_list)
        if len(png_list) == 0:
            return None

        out_mp4 = Path(out_mp4)
        frames_dir = out_mp4.parent / (out_mp4.stem + "_frames")
        if frames_dir.exists():
            shutil.rmtree(frames_dir)
        frames_dir.mkdir(parents=True, exist_ok=True)

        png_list = sorted(png_list)
        for i, p in enumerate(png_list):
            p = Path(p)
            link_path = frames_dir / f"{i:06d}.png"
            try:
                os.symlink(p, link_path)
            except Exception:
                shutil.copy2(p, link_path)

        cmd = [
            "ffmpeg", "-y",
            "-framerate", str(int(fps)),
            "-i", str(frames_dir / "%06d.png"),
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-movflags", "+faststart",
            str(out_mp4),
        ]
        subprocess.check_call(cmd)
        return str(out_mp4)

    def _make_video_imageio(png_list, out_gif, fps=8):
        import imageio.v2 as imageio
        imgs = [imageio.imread(p) for p in png_list]
        out_gif = str(out_gif)
        imageio.mimsave(out_gif, imgs, duration=1.0/max(int(fps), 1))
        return out_gif

    # -------------------------
    # paths
    # -------------------------
    step4_root = Path(step4_root).expanduser()
    if step5b_root is not None:
        step5b_root = Path(step5b_root).expanduser()

    if out_root is None:
        event_tag = _infer_event_tag_from_step4_root(step4_root)
        out_root = _default_out_root("step_iocentroid_outputs", event_tag)
    out_root = Path(out_root).expanduser()

    sbs = _find_sbs(step4_root)
    if len(sbs) == 0:
        raise ValueError(f"No SBxxx directories found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    # -------------------------
    # UI widgets
    # -------------------------
    source_dd = w.Dropdown(options=["Step4 raw", "Step5b corrected"], value="Step4 raw", description="Source")
    sb_dd = w.Dropdown(options=sbs, value=default_sb_ms, description="SB")
    refresh_btn = w.Button(description="Refresh FITS", button_style="")

    fits_sel = w.SelectMultiple(options=[], description="FITS", layout=w.Layout(width="95%", height="220px"))
    select_all_cb = w.Checkbox(value=False, description="Select ALL image.fits")
    overwrite_cb = w.Checkbox(value=False, description="Overwrite outputs")

    make_video_cb = w.Checkbox(value=False, description="Make video (mp4/gif)")
    fps_in = w.IntText(value=8, description="FPS")

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop hw (as)")
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low%")
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high%")
    contours_in = w.Text(value=",".join(map(str, default_contours)), description="Contours%")

    xmin_in = w.FloatText(value=float(default_roi[0]), description="xmin")
    xmax_in = w.FloatText(value=float(default_roi[1]), description="xmax")
    ymin_in = w.FloatText(value=float(default_roi[2]), description="ymin")
    ymax_in = w.FloatText(value=float(default_roi[3]), description="ymax")

    thresh_in = w.FloatText(value=float(default_thresh_frac), description="thresh_frac")
    minpts_in = w.IntText(value=int(default_min_points), description="min_points")

    show_inline_cb = w.Checkbox(value=True, description="Show images inline")
    max_inline_in = w.IntText(value=12, description="Max images")

    run_btn = w.Button(description="Run Step-5C (centroid)", button_style="success")
    out = w.Output()

    # -------------------------
    # refresh logic
    # -------------------------
    def _refresh(_=None):
        sb = sb_dd.value
        if source_dd.value == "Step4 raw":
            files = _list_step4_image_fits(step4_root, sb)
        else:
            files = _list_step5b_corr_fits(step5b_root, sb) if step5b_root is not None else []
        files = _load_sort_fits(files)
        fits_sel.options = files
        with out:
            clear_output()
            print("Source:", source_dd.value)
            print("SB:", sb)
            print("Step4 root:", step4_root)
            print("Step5b root:", step5b_root if step5b_root is not None else "(None)")
            print("Found:", len(files), "fits")
            if len(files) > 0:
                print("First:", os.path.basename(files[0]))
                print("Last :", os.path.basename(files[-1]))

    # -------------------------
    # run logic
    # -------------------------
    def _on_run(_):
        with out:
            clear_output()

            sb = sb_dd.value
            tag = _sb_tag(sb)
            source = source_dd.value

            # select files
            if select_all_cb.value:
                files = list(fits_sel.options)
            else:
                files = list(fits_sel.value)

            if len(files) == 0:
                print("No FITS selected.")
                return

            # params
            crop_hw = float(crop_in.value)
            clim = (float(clim_low.value), float(clim_high.value))
            contours = _parse_contours(contours_in.value)
            roi = (float(xmin_in.value), float(xmax_in.value), float(ymin_in.value), float(ymax_in.value))
            thresh_frac = float(thresh_in.value)
            min_points = int(minpts_in.value)

            sb_out = _ensure_dir(out_root / tag)
            qdir = _ensure_dir(sb_out / "quicklook_centroid")
            log_path = sb_out / "step5c_centroid.log"
            csv_path = sb_out / "centroid_results.csv"
            jsonl_path = sb_out / "centroid_results.jsonl"

            print("Source:", source)
            print("SB:", sb, "->", tag)
            print("Out:", sb_out)
            print("Quicklooks:", qdir)
            print("Overwrite:", overwrite_cb.value)
            print("Make video:", make_video_cb.value, "FPS:", int(fps_in.value))
            print("ROI (as):", roi)
            print("Fit params:", f"thresh_frac={thresh_frac}", f"min_points={min_points}")
            print("----")

            # write headers if new or overwrite
            if overwrite_cb.value or (not csv_path.exists()):
                with open(csv_path, "w") as f:
                    f.write("event_tag,sb,fname,tindex,obstime_isot,freq_mhz,roi_xmin,roi_xmax,roi_ymin,roi_ymax,thresh_frac,min_points,cen_Tx_as,cen_Ty_as,ok,png\n")
            if overwrite_cb.value and jsonl_path.exists():
                jsonl_path.unlink()

            pngs = []
            n_ok = 0

            event_tag = _infer_event_tag_from_step4_root(step4_root)

            for i, fpath in enumerate(files):
                fname = os.path.basename(fpath)
                tindex = _parse_tindex(fname)
                out_png = qdir / (fname.replace(".fits", "_centroid.png"))

                try:
                    png, meta = _quicklook_centroid_one(
                        fits_path=fpath,
                        out_png=str(out_png),
                        crop_hw_as=crop_hw,
                        clim_pct=clim,
                        contours=contours,
                        roi=roi,
                        thresh_frac=thresh_frac,
                        min_points=min_points,
                        overwrite=overwrite_cb.value,
                    )
                    pngs.append(png)
                    ok = bool(meta["ok"]) if meta else False
                    n_ok += 1 if ok else 0

                    # append CSV row
                    with open(csv_path, "a") as f:
                        f.write(
                            f"{event_tag},{tag},{fname},{tindex},{meta['obstime_isot']},{meta['freq_mhz']},"
                            f"{roi[0]},{roi[1]},{roi[2]},{roi[3]},"
                            f"{thresh_frac},{min_points},"
                            f"{meta['cen_Tx_as']},{meta['cen_Ty_as']},{int(ok)},{png}\n"
                        )

                    # append JSONL
                    rec = dict(
                        event_tag=event_tag,
                        sb=tag,
                        source=source,
                        fits=str(fpath),
                        fname=fname,
                        tindex=int(tindex),
                        obstime_isot=meta["obstime_isot"],
                        freq_mhz=meta["freq_mhz"],
                        roi=dict(xmin=roi[0], xmax=roi[1], ymin=roi[2], ymax=roi[3]),
                        thresh_frac=thresh_frac,
                        min_points=min_points,
                        cen_Tx_as=meta["cen_Tx_as"],
                        cen_Ty_as=meta["cen_Ty_as"],
                        ok=ok,
                        png=str(png),
                    )
                    with open(jsonl_path, "a") as f:
                        f.write(json.dumps(rec) + "\n")

                    # log
                    with open(log_path, "a") as f:
                        f.write(f"OK   {fpath} -> {png}\n")

                    # inline display control (avoid flooding)
                    if show_inline_cb.value and i < int(max_inline_in.value):
                        from IPython.display import Image as _Img, display as _disp
                        _disp(_Img(filename=png))

                except Exception as e:
                    with open(log_path, "a") as f:
                        f.write(f"FAIL {fpath} : {e}\n")
                    print("[FAIL]", fname, ":", e)

            print("Done:", f"{n_ok}/{len(files)} centroid fits OK")
            print("CSV :", csv_path)
            print("JSONL:", jsonl_path)
            print("Log :", log_path)

            # make video
            if make_video_cb.value and len(pngs) > 1:
                mp4 = sb_out / f"{tag}_centroid.mp4"
                gif = sb_out / f"{tag}_centroid.gif"
                try:
                    if shutil.which("ffmpeg"):
                        _make_video_ffmpeg(pngs, mp4, fps=int(fps_in.value))
                        print("Video:", mp4)
                    else:
                        raise RuntimeError("ffmpeg not found")
                except Exception:
                    try:
                        _make_video_imageio(pngs, gif, fps=int(fps_in.value))
                        print("GIF:", gif)
                    except Exception as e:
                        print("[WARN] video failed:", e)

    refresh_btn.on_click(_refresh)
    sb_dd.observe(_refresh, names="value")
    source_dd.observe(_refresh, names="value")
    run_btn.on_click(_on_run)

    _refresh()

    display(
        w.VBox([
            w.HBox([source_dd, sb_dd, refresh_btn, select_all_cb, overwrite_cb, make_video_cb, fps_in]),
            fits_sel,
            w.HBox([crop_in, clim_low, clim_high]),
            w.HBox([contours_in]),
            w.HBox([xmin_in, xmax_in, ymin_in, ymax_in]),
            w.HBox([thresh_in, minpts_in, show_inline_cb, max_inline_in]),
            run_btn,
            out
        ])
    )