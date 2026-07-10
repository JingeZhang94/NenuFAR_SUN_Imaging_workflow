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
      - Calibrator source selection
      - Calibrator event selection
      - Frequency range selection
      - Run scan_sun_and_cal_by_ymd

    New behavior:
      - supports AUTO / CAS_A / CYG_A / VIR_A
      - after successful run, writes latest_selected_plan.txt
    """
    import re
    from pathlib import Path

    import ipywidgets as w
    from IPython.display import display

    from nenufar_sb_scan import (
        scan_sun_and_cal_by_ymd,   # <-- new backend function
        list_event_dirs_by_date,
        pick_event_dir,
        pick_closest_calibrator,
    )

    BASE = Path(BASE)
    WORKROOT = Path(WORKROOT)

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

    cal_source = w.Dropdown(
        options=["AUTO", "CAS_A", "CYG_A", "VIR_A"],
        value="AUTO",
        description="Cal src",
    )

    cal = w.Dropdown(options=[], description="Cal")

    allfreq = w.Checkbox(value=False, description="All freq")
    fmin = w.FloatText(value=70.0, description="Fmin")
    fmax = w.FloatText(value=80.0, description="Fmax")

    preview = w.Dropdown(options=[30, 100, 300, "All"], value=30, description="Show")
    search = w.Text(value="", description="Search", placeholder="e.g. SB407 / 407 / 60.3")

    btn = w.Button(description="Run", button_style="success")
    out = w.Output()

    # ---------- helpers ----------
    def _apply_search(df, q: str):
        q = (q or "").strip()
        if not q:
            return df
        qlow = q.lower()
        mask = df.astype(str).apply(lambda row: row.str.lower().str.contains(qlow, na=False)).any(axis=1)
        return df.loc[mask]

    def _find_latest_plan(workroot: Path):
        cands = sorted(
            workroot.glob("*/selected_sb_pair_list.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        return cands[0] if cands else None

    def _write_latest_plan_file(workroot: Path, plan_path: Path):
        latest_txt = workroot / "latest_selected_plan.txt"
        latest_txt.write_text(str(plan_path) + "\n")
        return latest_txt

    def _get_calibrator_candidates_for_date(y, m, d, source_name):
        """
        Return candidate calibrator event dirs for a given date and source.
        source_name:
          - AUTO
          - CAS_A / CYG_A / VIR_A
        """
        if source_name == "AUTO":
            all_cands = []
            for tag in ["CAS_A_TRACKING", "CYG_A_TRACKING", "VIR_A_TRACKING"]:
                all_cands.extend(list_event_dirs_by_date(str(BASE), y, m, d, tag))
            return sorted(all_cands)

        track_tag = f"{source_name}_TRACKING"
        return list_event_dirs_by_date(str(BASE), y, m, d, track_tag)

    # ---------- refresh calibrator ----------
    def refresh_cal(*_):
        try:
            sun_ev = pick_event_dir(str(BASE), year.value, month.value, day.value, "SUN_TRACKING")
            cands = _get_calibrator_candidates_for_date(
                year.value, month.value, day.value, cal_source.value
            )
            cal.options = cands
            if cands:
                cal.value = pick_closest_calibrator(sun_ev, cands)
            else:
                cal.value = None
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

    def on_cal_source_change(change):
        if change.get("name") != "value":
            return
        refresh_cal()

    year.observe(on_year_change, names="value")
    month.observe(on_month_change, names="value")
    day.observe(on_day_change, names="value")
    cal_source.observe(on_cal_source_change, names="value")

    # ---------- run ----------
    def run(_):
        with out:
            out.clear_output()

            fr = None if allfreq.value else (float(fmin.value), float(fmax.value))

            df_sel, meta = scan_sun_and_cal_by_ymd(
                str(BASE),
                str(WORKROOT),
                year.value,
                month.value,
                day.value,
                fr,
                cal_source=cal_source.value,
                cal_event_dir=cal.value,
            )

            df_show = _apply_search(df_sel, search.value)

            print(meta)
            print(f"Found {len(df_sel)} SBs total; after search filter: {len(df_show)}")
            print(f"Showing: {preview.value}")

            latest_plan = _find_latest_plan(WORKROOT)
            if latest_plan is not None:
                latest_txt = _write_latest_plan_file(WORKROOT, latest_plan)
                print(f"\n[WRITE] latest selected plan -> {latest_txt}")
                print(f"[PLAN ] {latest_plan}")
            else:
                print("\n[WARN] No selected_sb_pair_list.json found under WORKROOT after run.")

            if preview.value == "All":
                display(df_show)
            else:
                display(df_show.head(int(preview.value)))

    btn.on_click(run)

    display(
        w.HBox([year, month, day, cal_source, cal]),
        w.HBox([allfreq, fmin, fmax, preview, search, btn]),
        out,
    )
def run_step1_ui(plan_json_path, out_root=None):
    """
    Small UI to run DP3 Step-1 (prepare Cal + Sun) for selected SBs.

    Inputs
    - plan_json_path: selected_sb_pair_list.json
        must include: selected_sb (or selected_sbs) + sun_ms
        supports both old casa_* fields and new cal_* fields
    - out_root: where to write outputs; default: alongside plan json folder

    What it runs (per SB)
    - Calibrator branch(es):
        1A) AOFlagger + averager
        1B) preflagger (in-place)
    - Solar branch:
        1C) clear flags

    Modes
    - dropdown : use chosen calibrator from JSON if present
    - pre      : prepare only pre-calibrator + sun
    - post     : prepare only post-calibrator + sun
    - both     : prepare both pre- and post-calibrators + sun

    New behavior
    - supports generic calibrator source names from JSON:
        CAL_SOURCE_CHOSEN = CAS_A / CYG_A / VIR_A
    - output file names follow calibrator source:
        CasApre_SB359_prep.MS / CygApre_SB359_prep.MS / VirApost_SB359_prep.MS
    - if pre/post is missing, it skips gracefully
    """
    import json
    import re
    import shutil
    import shlex
    import subprocess
    from pathlib import Path

    import ipywidgets as w
    from IPython.display import display

    # ---- container command ----
    SIF = "/home/jzhang/LOFARimgCode/linc_latest.sif"
    BIND = "-B /databf -B /data"
    ENGINE = "apptainer"
    DP3_CMD = f"{ENGINE} exec {BIND} {SIF} DP3"

    plan_json_path = Path(plan_json_path)
    payload = json.load(open(plan_json_path, "r"))

    # --------------------------------------------------
    # payload core arrays
    # --------------------------------------------------
    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    sun_ms_list = payload.get("sun_ms") or []

    # new generic fields preferred
    cal_ms_list = payload.get("cal_chosen_ms") or []
    cal_pre_ms_list = payload.get("cal_pre_ms") or []
    cal_post_ms_list = payload.get("cal_post_ms") or []

    # backward compatibility
    if not cal_ms_list:
        cal_ms_list = payload.get("casa_ms") or []
    if not cal_pre_ms_list:
        cal_pre_ms_list = payload.get("casa_pre_ms") or payload.get("casa_pre") or []
    if not cal_post_ms_list:
        cal_post_ms_list = payload.get("casa_post_ms") or payload.get("casa_post") or []

    cal_source_chosen = payload.get("cal_source_chosen") or "CAS_A"

    if not (sbs and sun_ms_list):
        raise ValueError("JSON must contain 'selected_sb' (or 'selected_sbs') and 'sun_ms'.")

    # --------------------------------------------------
    # helpers for calibrator naming
    # --------------------------------------------------
    def _cal_display_tag(cal_source):
        mapping = {
            "CAS_A": "CasA",
            "CYG_A": "CygA",
            "VIR_A": "VirA",
        }
        return mapping.get(str(cal_source).upper(), str(cal_source))

    CAL_TAG = _cal_display_tag(cal_source_chosen)

    # mappings
    sun_map = dict(zip(sbs, sun_ms_list))
    cal_map = dict(zip(sbs, cal_ms_list)) if cal_ms_list else {}
    cal_pre_map = dict(zip(sbs, cal_pre_ms_list)) if cal_pre_ms_list else {}
    cal_post_map = dict(zip(sbs, cal_post_ms_list)) if cal_post_ms_list else {}

    # --------------------------------------------------
    # frequency helpers
    # --------------------------------------------------
    def _normalize_sb_key(sb):
        m = re.search(r"(SB\d+)", str(sb))
        return m.group(1) if m else str(sb)

    def _build_freq_map_from_payload(payload, sbs):
        sb_norm_to_full = {_normalize_sb_key(sb): sb for sb in sbs}
        freq_map = {}

        candidate_array_keys = [
            "ctr_mhz",
            "selected_freq_mhz", "selected_freqs_mhz",
            "freq_mhz", "freqs_mhz",
            "selected_freq", "selected_freqs",
            "frequency_mhz", "frequencies_mhz",
            "frequency", "frequencies",
        ]
        for key in candidate_array_keys:
            arr = payload.get(key)
            if isinstance(arr, list) and len(arr) == len(sbs):
                ok = True
                for sb, f in zip(sbs, arr):
                    try:
                        freq_map[sb] = float(f)
                    except Exception:
                        ok = False
                        break
                if ok and len(freq_map) == len(sbs):
                    return freq_map

        candidate_dict_keys = [
            "sb_freq_map", "sb_to_freq_mhz", "sb_frequency_map",
            "freq_map", "frequency_map",
        ]
        for key in candidate_dict_keys:
            d = payload.get(key)
            if isinstance(d, dict):
                for k, v in d.items():
                    sb_norm = _normalize_sb_key(k)
                    if sb_norm in sb_norm_to_full:
                        try:
                            freq_map[sb_norm_to_full[sb_norm]] = float(v)
                        except Exception:
                            pass
                if freq_map:
                    return freq_map

        candidate_record_keys = [
            "selected_sb_info", "selected_sb_records", "sb_info", "sb_records"
        ]
        for key in candidate_record_keys:
            recs = payload.get(key)
            if isinstance(recs, list):
                for rec in recs:
                    if not isinstance(rec, dict):
                        continue
                    sb_val = rec.get("sb") or rec.get("SB") or rec.get("selected_sb")
                    f_val = (
                        rec.get("freq_mhz") or rec.get("frequency_mhz") or
                        rec.get("freq") or rec.get("frequency")
                    )
                    if sb_val is None or f_val is None:
                        continue
                    sb_norm = _normalize_sb_key(sb_val)
                    if sb_norm in sb_norm_to_full:
                        try:
                            freq_map[sb_norm_to_full[sb_norm]] = float(f_val)
                        except Exception:
                            pass
                if freq_map:
                    return freq_map

        raw_selected = payload.get("selected_sb") or payload.get("selected_sbs")
        if isinstance(raw_selected, list) and raw_selected and isinstance(raw_selected[0], dict):
            for rec in raw_selected:
                sb_val = rec.get("sb") or rec.get("SB") or rec.get("selected_sb")
                f_val = (
                    rec.get("freq_mhz") or rec.get("frequency_mhz") or
                    rec.get("freq") or rec.get("frequency")
                )
                if sb_val is None or f_val is None:
                    continue
                sb_norm = _normalize_sb_key(sb_val)
                if sb_norm in sb_norm_to_full:
                    try:
                        freq_map[sb_norm_to_full[sb_norm]] = float(f_val)
                    except Exception:
                        pass
            if freq_map:
                return freq_map

        return freq_map

    freq_map = _build_freq_map_from_payload(payload, sbs)

    def _sb_label(sb):
        if sb in freq_map:
            return f"{sb}   |   {freq_map[sb]:.1f} MHz"
        return sb

    sb_options = [(_sb_label(sb), sb) for sb in sbs]

    # --------------------------------------------------
    # cal mode options
    # --------------------------------------------------
    cal_mode_opts = ["dropdown"]
    if cal_pre_ms_list and cal_post_ms_list:
        cal_mode_opts = ["dropdown", "pre", "post", "both"]
    elif cal_pre_ms_list:
        cal_mode_opts = ["dropdown", "pre"]
    elif cal_post_ms_list:
        cal_mode_opts = ["dropdown", "post"]

    cal_mode = w.ToggleButtons(options=cal_mode_opts, value="dropdown", description="Cal mode")

    sb_select = w.SelectMultiple(
        options=sb_options,
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

    # --------------------------------------------------
    # filesystem helpers
    # --------------------------------------------------
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

    def _normalize_ms_entry(entry, sb):
        """
        Accept either:
        - full MS path
        - event dir -> append /L1/SBxxx.MS
        - L1 dir    -> append /SBxxx.MS
        """
        if entry is None:
            return None
        p = Path(str(entry))
        if str(p).endswith(".MS") or p.name.endswith(".MS"):
            return str(p)
        if p.name == "L1":
            return str(p / sb)
        return str(p / "L1" / sb)

    def _get_dropdown_cal_ms_for_sb(sb):
        if not cal_map:
            return None
        return _normalize_ms_entry(cal_map.get(sb), sb)

    def _get_pre_cal_ms_for_sb(sb):
        if not cal_pre_map:
            return None
        return _normalize_ms_entry(cal_pre_map.get(sb), sb)

    def _get_post_cal_ms_for_sb(sb):
        if not cal_post_map:
            return None
        return _normalize_ms_entry(cal_post_map.get(sb), sb)

    # --------------------------------------------------
    # build jobs
    # --------------------------------------------------
    def _build_cal_jobs_for_sb(sb):
        """
        Return a list of calibrator jobs for this SB.
        Each job:
          - label: dropdown / pre / post
          - input_ms
          - output_ms
          - log_a
          - log_b
        """
        sb_tag = sb.replace(".MS", "")
        sb_dir = out_root / sb_tag
        jobs = []

        if cal_mode.value == "dropdown":
            cal_ms = _get_dropdown_cal_ms_for_sb(sb)
            if cal_ms:
                jobs.append(dict(
                    label="dropdown",
                    input_ms=cal_ms,
                    output_ms=sb_dir / f"{CAL_TAG}_{sb_tag}_prep.MS",
                    log_a=sb_dir / f"dp3_{CAL_TAG.lower()}_aoflag_avg.log",
                    log_b=sb_dir / f"dp3_{CAL_TAG.lower()}_preflag.log",
                ))

        elif cal_mode.value == "pre":
            cal_ms = _get_pre_cal_ms_for_sb(sb)
            if cal_ms:
                jobs.append(dict(
                    label="pre",
                    input_ms=cal_ms,
                    output_ms=sb_dir / f"{CAL_TAG}pre_{sb_tag}_prep.MS",
                    log_a=sb_dir / f"dp3_{CAL_TAG.lower()}_pre_aoflag_avg.log",
                    log_b=sb_dir / f"dp3_{CAL_TAG.lower()}_pre_preflag.log",
                ))

        elif cal_mode.value == "post":
            cal_ms = _get_post_cal_ms_for_sb(sb)
            if cal_ms:
                jobs.append(dict(
                    label="post",
                    input_ms=cal_ms,
                    output_ms=sb_dir / f"{CAL_TAG}post_{sb_tag}_prep.MS",
                    log_a=sb_dir / f"dp3_{CAL_TAG.lower()}_post_aoflag_avg.log",
                    log_b=sb_dir / f"dp3_{CAL_TAG.lower()}_post_preflag.log",
                ))

        elif cal_mode.value == "both":
            cal_pre_ms = _get_pre_cal_ms_for_sb(sb)
            cal_post_ms = _get_post_cal_ms_for_sb(sb)

            if cal_pre_ms:
                jobs.append(dict(
                    label="pre",
                    input_ms=cal_pre_ms,
                    output_ms=sb_dir / f"{CAL_TAG}pre_{sb_tag}_prep.MS",
                    log_a=sb_dir / f"dp3_{CAL_TAG.lower()}_pre_aoflag_avg.log",
                    log_b=sb_dir / f"dp3_{CAL_TAG.lower()}_pre_preflag.log",
                ))

            if cal_post_ms:
                jobs.append(dict(
                    label="post",
                    input_ms=cal_post_ms,
                    output_ms=sb_dir / f"{CAL_TAG}post_{sb_tag}_prep.MS",
                    log_a=sb_dir / f"dp3_{CAL_TAG.lower()}_post_aoflag_avg.log",
                    log_b=sb_dir / f"dp3_{CAL_TAG.lower()}_post_preflag.log",
                ))

        return jobs

    # --------------------------------------------------
    # run
    # --------------------------------------------------
    def on_run(_):
        with out:
            out.clear_output()

            chosen = list(sb_select.value)
            if not chosen:
                print("No SB selected.")
                return

            print("Plan:", plan_json_path)
            print("Out root:", out_root)
            print("Cal source chosen:", cal_source_chosen)
            print("Cal tag:", CAL_TAG)
            print("Cal mode:", cal_mode.value)
            print("Overwrite:", overwrite.value)
            print("SBs:", ", ".join([_sb_label(sb) for sb in chosen]))

            for sb in chosen:
                if sb not in sun_map:
                    print(f"[WARN] SB not in sun_map: {sb} (skip)")
                    continue

                sun_ms = sun_map[sb]
                sb_tag = sb.replace(".MS", "")
                sb_dir = out_root / sb_tag
                sb_dir.mkdir(parents=True, exist_ok=True)

                cal_jobs = _build_cal_jobs_for_sb(sb)

                if len(cal_jobs) == 0:
                    print(f"[WARN] No calibrator MS found for {sb} under cal_mode={cal_mode.value} (skip SB)")
                    continue

                sun_out = sb_dir / f"SUN_{sb_tag}_prep.MS"
                log_sun = sb_dir / "dp3_sun_clearflags.log"

                if overwrite.value:
                    for job in cal_jobs:
                        _rm_tree(job["output_ms"])
                    _rm_tree(sun_out)

                # ---------------- Calibrator branches ----------------
                for job in cal_jobs:
                    cal_in = job["input_ms"]
                    cal_out = job["output_ms"]
                    log1 = job["log_a"]
                    log2 = job["log_b"]
                    label = job["label"]

                    cmd1 = (
                        f"{DP3_CMD} msin={shlex.quote(cal_in)} "
                        f"msout={shlex.quote(str(cal_out))} "
                        f"steps=[flag,averager] "
                        f"flag.type=aoflagger "
                        f"averager.timestep=1 averager.freqstep=1"
                    )

                    cmd2 = (
                        f"{DP3_CMD} msin={shlex.quote(str(cal_out))} "
                        f"msout=. "
                        f"steps=[flag] "
                        f"flag.type=preflagger "
                        f"flag.baseline='MR100NEN&&*;MR101NEN&&*;MR102NEN&&*;MR103NEN&&*'"
                    )

                    print(f"\n=== Step-1A (Cal {label}): {_sb_label(sb)} ===")
                    if (not overwrite.value) and cal_out.exists():
                        print(f"[SKIP 1A {label}] exists: {cal_out}")
                    else:
                        print(cmd1)
                        _run_to_log(cmd1, log1)

                    print(f"=== Step-1B (Cal {label} preflag): {_sb_label(sb)} ===")
                    if not cal_out.exists():
                        raise RuntimeError(f"1B ({label}) needs {cal_out}, but it does not exist.")
                    if (not overwrite.value) and log2.exists():
                        print(f"[SKIP 1B {label}] log exists: {log2}")
                    else:
                        print(cmd2)
                        _run_to_log(cmd2, log2)

                # ---------------- Solar branch ----------------
                cmd3 = (
                    f"{DP3_CMD} msin={shlex.quote(sun_ms)} "
                    f"msout={shlex.quote(str(sun_out))} "
                    f"steps=[flag] "
                    f"flag.type=preflagger flag.mode=clear flag.baseline='*&&*'"
                )

                print(f"=== Step-1C (Sun prep): {_sb_label(sb)} ===")
                if (not overwrite.value) and sun_out.exists():
                    print(f"[SKIP 1C] exists: {sun_out}")
                else:
                    print(cmd3)
                    _run_to_log(cmd3, log_sun)

                print(f"OK: {_sb_label(sb)} -> {sb_dir}")

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
                  SUN_{SB}_prep.MS exists, only those available SBs will be shown in the UI,
                  and these prepared SUN MS will be used as msin.
    - out_root: where to write step-2 outputs. default: sibling folder of plan
    """

    import json
    import shutil
    from pathlib import Path
    import ipywidgets as w
    from IPython.display import display
    import subprocess
    import shlex

    # ---- container / dp3 ----
    SIF = "/home/jzhang/LOFARimgCode/linc_latest.sif"
    BIND = "-B /databf -B /data"
    ENGINE = "apptainer"
    DP3_CMD = f"{ENGINE} exec {BIND} {SIF} DP3"

    plan_json_path = Path(plan_json_path)
    payload = json.load(open(plan_json_path, "r"))

    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    sun_ms_list = payload.get("sun_ms") or []
    ctr_mhz = payload.get("ctr_mhz") or []

    if not (sbs and sun_ms_list):
        raise ValueError("JSON must contain 'selected_sb' (or selected_sbs) and 'sun_ms'.")

    sun_map = dict(zip(sbs, sun_ms_list))
    freq_map = dict(zip(sbs, ctr_mhz)) if len(ctr_mhz) == len(sbs) else {}

    # defaults
    if out_root is None:
        out_root = plan_json_path.parent.parent / "step2_outputs"
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    step1_root = Path(step1_root) if step1_root else None

    # --------------------------------------------------
    # helpers
    # --------------------------------------------------
    def _sb_label(sb):
        if sb in freq_map:
            try:
                return f"{sb}   |   {float(freq_map[sb]):.1f} MHz"
            except Exception:
                return sb
        return sb

    def _step1_sun_prep_path(sb: str):
        if step1_root is None:
            return None
        sb_tag = sb.replace(".MS", "")
        return step1_root / sb_tag / f"SUN_{sb_tag}_prep.MS"

    def _available_sbs():
        """
        If step1_root is provided, only show SBs whose Step-1 prepared SUN MS exists.
        Otherwise, fall back to all SBs from the plan.
        """
        if step1_root is None:
            return list(sbs)

        good = []
        for sb in sbs:
            cand = _step1_sun_prep_path(sb)
            if cand is not None and cand.exists():
                good.append(sb)
        return good

    available_sbs = _available_sbs()

    if len(available_sbs) == 0:
        raise FileNotFoundError(
            f"No available SUN_*_prep.MS found under step1_root={step1_root}. "
            f"Please run Step-1 first, or check the step1_root path."
        )

    sb_options = [(_sb_label(sb), sb) for sb in available_sbs]

    # ---- widgets ----
    default_value = ("SB359.MS",) if "SB359.MS" in available_sbs else (available_sbs[0],)

    sb_select = w.SelectMultiple(
        options=sb_options,
        value=default_value,
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
        cand = _step1_sun_prep_path(sb)
        if cand is not None and cand.exists():
            return str(cand)
        return str(sun_map[sb])

    def _safe_remove_ms(path: Path):
        if path.exists():
            shutil.rmtree(path)

    # --------------------------------------------------
    # run
    # --------------------------------------------------
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
            print(f"Available SBs from Step-1: {len(available_sbs)}")
            print("SBs:", ", ".join([_sb_label(sb) for sb in chosen]))

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
                    print(f"\n[Skip] ROI exists for {_sb_label(sb)}: {msout}")
                    continue

                cmd_roi = (
                    f"{DP3_CMD} "
                    f"msin={shlex.quote(str(msin))} "
                    f"msout={shlex.quote(str(msout))} "
                    f"steps=[] "
                    f"msin.starttime={shlex.quote(t_start.value)} "
                    f"msin.endtime={shlex.quote(t_end.value)}"
                )

                print(f"\n=== ROI cut {_sb_label(sb)} ===")
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
    Step-3: Calibration

    Upgraded version:
    - show SB + frequency in UI
    - only show SBs that already have Step-2 ROI products
    - support generic calibrator source from JSON:
        CAS_A / CYG_A / VIR_A
    - support skymodel path passed via `sourcedb` argument
    - resolve Step-1 calibrator prep MS using new naming:
        CygApre_SB359_prep.MS / VirApost_SB359_prep.MS / CasA_SB359_prep.MS
      while keeping backward compatibility with old CasA-only names
    - IMPORTANT FIX:
        single-cal mode now shifts the instrument STARTY/ENDY to the Sun ROI time window
        before applycal, so the solution really applies to Sun data.
    - both mode:
        * run gaincal on pre anchor
        * run gaincal on post anchor
        * stop before interpolation/applycal (for now)

    Modes
    - dropdown : single-calibrator path -> {CalTag}_{SB}_prep.MS
    - pre      : single-calibrator path -> {CalTag}pre_{SB}_prep.MS
    - post     : single-calibrator path -> {CalTag}post_{SB}_prep.MS
    - both     : dual-anchor mode, run both gaincals only
    """

    import json
    import time
    import shlex
    import shutil
    import subprocess
    from pathlib import Path

    import ipywidgets as w
    from IPython.display import display

    # --------------------------------------------------
    # paths
    # --------------------------------------------------
    plan_json_path = Path(plan_json_path)
    step1_root = Path(step1_root)
    step2_root = Path(step2_root)
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    sourcedb = Path(sourcedb)
    if not sourcedb.exists():
        raise FileNotFoundError(f"Calibration model path not found: {sourcedb}")

    DP3_CMD = f"{engine} exec {bind} {shlex.quote(sif)} DP3"

    payload = json.load(open(plan_json_path, "r"))

    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    ctr_mhz = payload.get("ctr_mhz") or []

    if not sbs:
        raise ValueError("JSON must contain 'selected_sb' or 'selected_sbs'.")

    freq_map = dict(zip(sbs, ctr_mhz)) if len(ctr_mhz) == len(sbs) else {}

    # new generic fields preferred
    cal_pre_list = payload.get("cal_pre_ms") or []
    cal_post_list = payload.get("cal_post_ms") or []
    cal_ms_list = payload.get("cal_chosen_ms") or []

    # backward compatibility
    if not cal_pre_list:
        cal_pre_list = payload.get("casa_pre_ms") or payload.get("casa_pre") or []
    if not cal_post_list:
        cal_post_list = payload.get("casa_post_ms") or payload.get("casa_post") or []
    if not cal_ms_list:
        cal_ms_list = payload.get("casa_ms") or []

    cal_source_chosen = payload.get("cal_source_chosen") or "CAS_A"

    # --------------------------------------------------
    # source/tag helpers
    # --------------------------------------------------
    def _cal_display_tag(cal_source):
        mapping = {
            "CAS_A": "CasA",
            "CYG_A": "CygA",
            "VIR_A": "VirA",
        }
        return mapping.get(str(cal_source).upper(), str(cal_source))

    def _cal_gaincal_source_name(cal_source):
        mapping = {
            "CAS_A": "CasA",
            "CYG_A": "CygA",
            "VIR_A": "VirA",
        }
        return mapping.get(str(cal_source).upper(), str(cal_source))

    CAL_TAG = _cal_display_tag(cal_source_chosen)
    GAINCAL_SOURCE_NAME = _cal_gaincal_source_name(cal_source_chosen)

    # --------------------------------------------------
    # helper: display labels
    # --------------------------------------------------
    def _sb_label(sb):
        if sb in freq_map:
            try:
                return f"{sb}   |   {float(freq_map[sb]):.1f} MHz"
            except Exception:
                return sb
        return sb

    # --------------------------------------------------
    # helper: available SBs from Step-2
    # --------------------------------------------------
    def _step2_roi_path(sb: str):
        sb_tag = sb.replace(".MS", "")
        return step2_root / "ROI" / sb_tag / f"SUN_{sb_tag}_ROI.MS"

    def _available_sbs():
        good = []
        for sb in sbs:
            if _step2_roi_path(sb).exists():
                good.append(sb)
        return good

    available_sbs = _available_sbs()
    if len(available_sbs) == 0:
        raise FileNotFoundError(
            f"No available SUN_*_ROI.MS found under step2_root={step2_root}. "
            f"Please run Step-2 first, or check the step2_root path."
        )

    sb_options = [(_sb_label(sb), sb) for sb in available_sbs]

    # --------------------------------------------------
    # cal mode options
    # --------------------------------------------------
    cal_mode_opts = ["dropdown"]
    if cal_pre_list and cal_post_list:
        cal_mode_opts = ["dropdown", "pre", "post", "both"]
    elif cal_pre_list:
        cal_mode_opts = ["dropdown", "pre"]
    elif cal_post_list:
        cal_mode_opts = ["dropdown", "post"]

    cal_mode = w.ToggleButtons(options=cal_mode_opts, value="dropdown", description="Cal mode")
    overwrite = w.Checkbox(value=False, description="Overwrite existing outputs")

    default_value = ("SB359.MS",) if "SB359.MS" in available_sbs else (available_sbs[0],)

    sb_select = w.SelectMultiple(
        options=sb_options,
        value=default_value,
        description="SBs",
        layout=w.Layout(width="95%", height="180px"),
    )

    run_btn = w.Button(description="Run Step-3 (Calib)", button_style="success")
    out = w.Output()

    # --------------------------------------------------
    # low-level helpers
    # --------------------------------------------------
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

    def _find_existing_ms(candidates):
        for p in candidates:
            p = Path(p)
            if p.exists():
                return p
        return None

    def _copy_and_shift_instrument_to_target_ms(inst_in: Path, target_ms: Path, inst_out: Path, pad_s: float = 1.0, overwrite: bool = False):
        """
        Copy an instrument table and shift all STARTY/ENDY rows to cover the TIME range of target_ms.
        """
        import numpy as np
        from casacore.tables import table

        inst_in = Path(inst_in)
        target_ms = Path(target_ms)
        inst_out = Path(inst_out)

        if not inst_in.exists():
            raise FileNotFoundError(f"Input instrument not found: {inst_in}")
        if not target_ms.exists():
            raise FileNotFoundError(f"Target MS not found: {target_ms}")

        if inst_out.exists():
            if overwrite:
                shutil.rmtree(inst_out)
            else:
                return inst_out

        shutil.copytree(inst_in, inst_out)

        # read target TIME range
        t_ms = table(str(target_ms))
        times = t_ms.getcol("TIME")
        t_ms.close()

        tmin = float(np.min(times))
        tmax = float(np.max(times))

        # write shifted STARTY / ENDY
        t_inst = table(str(inst_out), readonly=False)
        starty = t_inst.getcol("STARTY")
        endy = t_inst.getcol("ENDY")

        starty[:] = tmin - pad_s
        endy[:] = tmax + pad_s

        t_inst.putcol("STARTY", starty)
        t_inst.putcol("ENDY", endy)
        t_inst.flush()
        t_inst.close()

        return inst_out

    # --------------------------------------------------
    # step1 product resolvers
    # --------------------------------------------------
    def _resolve_calib_ms_single(sb: str):
        sb_tag = sb.replace(".MS", "")
        sb_dir = step1_root / sb_tag

        if cal_mode.value == "pre":
            candidates = [
                sb_dir / f"{CAL_TAG}pre_{sb_tag}_prep.MS",
                sb_dir / f"{CAL_TAG}_{sb_tag}_prep.MS",
                sb_dir / f"CasApre_{sb_tag}_prep.MS",   # backward compatibility
            ]
            return _find_existing_ms(candidates)

        elif cal_mode.value == "post":
            candidates = [
                sb_dir / f"{CAL_TAG}post_{sb_tag}_prep.MS",
                sb_dir / f"{CAL_TAG}_{sb_tag}_prep.MS",
                sb_dir / f"CasApost_{sb_tag}_prep.MS",  # backward compatibility
            ]
            return _find_existing_ms(candidates)

        else:  # dropdown
            candidates = [
                sb_dir / f"{CAL_TAG}_{sb_tag}_prep.MS",
                sb_dir / f"CasA_{sb_tag}_prep.MS",      # backward compatibility
            ]
            return _find_existing_ms(candidates)

    def _resolve_calib_ms_both(sb: str):
        sb_tag = sb.replace(".MS", "")
        sb_dir = step1_root / sb_tag

        pre_candidates = [
            sb_dir / f"{CAL_TAG}pre_{sb_tag}_prep.MS",
            sb_dir / f"{CAL_TAG}_{sb_tag}_prep.MS",
            sb_dir / f"CasApre_{sb_tag}_prep.MS",
        ]
        post_candidates = [
            sb_dir / f"{CAL_TAG}post_{sb_tag}_prep.MS",
            sb_dir / f"{CAL_TAG}_{sb_tag}_prep.MS",
            sb_dir / f"CasApost_{sb_tag}_prep.MS",
        ]

        return _find_existing_ms(pre_candidates), _find_existing_ms(post_candidates)

    # --------------------------------------------------
    # parset builders
    # --------------------------------------------------
    def _gaincal_parset_text(calib_ms: Path):
        return f"""
msin={calib_ms}
msout=.

steps=[gaincal]

gaincal.usebeammodel=False

# Solve one gain over the full calibrator time range
gaincal.solint=0
# Solve one gain over the full channel range
gaincal.nchan=0

gaincal.sources={GAINCAL_SOURCE_NAME}
gaincal.sourcedb={sourcedb}
gaincal.onebeamperpatch=False
gaincal.caltype=diagonal
""".lstrip()

    def _applycal_parset_text(sun_roi_ms: Path, parmdb: Path):
        return f"""
msin={sun_roi_ms}
msout=.
msin.datacolumn=DATA
msout.datacolumn=CORR_NO_BEAM

steps=[applycal]

applycal.parmdb={parmdb}
applycal.correction=gain
applycal.updateweights=True
""".lstrip()

    # --------------------------------------------------
    # run
    # --------------------------------------------------
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
            print("Cal source chosen:", cal_source_chosen)
            print("Cal tag:", CAL_TAG)
            print("Gaincal source name:", GAINCAL_SOURCE_NAME)
            print("Cal mode:", cal_mode.value)
            print("Overwrite:", overwrite.value)
            print(f"Available SBs from Step-2: {len(available_sbs)}")
            print("SBs:", ", ".join([_sb_label(sb) for sb in chosen]))
            print("Using calibration model:", sourcedb)
            print()

            t0 = time.time()

            for sb in chosen:
                sb_tag = sb.replace(".MS", "")
                sb_dir = out_root / sb_tag
                sb_dir.mkdir(parents=True, exist_ok=True)

                # Step-2 ROI product
                sun_roi_ms = _step2_roi_path(sb)
                if not sun_roi_ms.exists():
                    raise FileNotFoundError(f"Missing Step-2 Sun ROI MS for {sb}: {sun_roi_ms}")

                # ==========================================================
                # both mode: dual gain anchors only (keep for later interpolation)
                # ==========================================================
                if cal_mode.value == "both":
                    calib_pre_ms, calib_post_ms = _resolve_calib_ms_both(sb)

                    if calib_pre_ms is None:
                        raise FileNotFoundError(f"Missing Step-1 pre calibrator MS for {sb}")
                    if calib_post_ms is None:
                        raise FileNotFoundError(f"Missing Step-1 post calibrator MS for {sb}")

                    parmdb_pre = calib_pre_ms / "instrument"
                    parmdb_post = calib_post_ms / "instrument"

                    par_gain_pre = sb_dir / f"{sb_tag}_gaincal_pre.parset"
                    par_gain_post = sb_dir / f"{sb_tag}_gaincal_post.parset"

                    gaincal_pre_text = _gaincal_parset_text(calib_pre_ms)
                    gaincal_post_text = _gaincal_parset_text(calib_post_ms)

                    wrote_pre = _write_text(par_gain_pre, gaincal_pre_text, force=overwrite.value)
                    wrote_post = _write_text(par_gain_post, gaincal_post_text, force=overwrite.value)

                    if wrote_pre:
                        print(f"[WRITE parset] {par_gain_pre}")
                    if wrote_post:
                        print(f"[WRITE parset] {par_gain_post}")

                    log_gain_pre = sb_dir / "01_gaincal_pre.log"
                    log_gain_post = sb_dir / "02_gaincal_post.log"

                    if parmdb_pre.exists() and (not overwrite.value):
                        print(f"[SKIP gaincal pre] {_sb_label(sb)}: parmdb exists ({parmdb_pre})")
                    else:
                        print(f"\n=== Step-3A gaincal PRE: {_sb_label(sb)} ===")
                        cmd_gain_pre = f"{DP3_CMD} {shlex.quote(str(par_gain_pre))}"
                        print(cmd_gain_pre)
                        _run_to_log(cmd_gain_pre, log_gain_pre)

                    if parmdb_post.exists() and (not overwrite.value):
                        print(f"[SKIP gaincal post] {_sb_label(sb)}: parmdb exists ({parmdb_post})")
                    else:
                        print(f"=== Step-3B gaincal POST: {_sb_label(sb)} ===")
                        cmd_gain_post = f"{DP3_CMD} {shlex.quote(str(par_gain_post))}"
                        print(cmd_gain_post)
                        _run_to_log(cmd_gain_post, log_gain_post)

                    print(f"[INFO] Both-mode currently stops after generating PRE/POST gain anchors.")
                    print(f"[INFO] Interpolated applycal is not implemented yet for {_sb_label(sb)}.")
                    print(f"OK: {_sb_label(sb)} -> {sb_dir}\n")
                    continue

                # ==========================================================
                # single-calibrator mode
                # FIXED: shift instrument time window to Sun ROI before applycal
                # ==========================================================
                calib_ms = _resolve_calib_ms_single(sb)
                if calib_ms is None:
                    raise FileNotFoundError(
                        f"Missing Step-1 calibrator MS for {sb}. "
                        f"Tried CAL_TAG={CAL_TAG} under mode={cal_mode.value}"
                    )

                parmdb_raw = calib_ms / "instrument"
                parmdb_shifted = sb_dir / "instrument_sunwindow"

                par_gain = sb_dir / f"{sb_tag}_gaincal.parset"
                par_apply = sb_dir / f"{sb_tag}_applycal.parset"

                gaincal_parset = _gaincal_parset_text(calib_ms)

                # Step-3A: gaincal on calibrator prep MS
                wrote_gain = _write_text(par_gain, gaincal_parset, force=overwrite.value)
                if wrote_gain:
                    print(f"[WRITE parset] {par_gain}")

                log_gain = sb_dir / "01_gaincal.log"
                if parmdb_raw.exists() and (not overwrite.value):
                    print(f"[SKIP gaincal] {_sb_label(sb)}: parmdb exists ({parmdb_raw})")
                else:
                    print(f"\n=== Step-3A gaincal: {_sb_label(sb)} ({cal_mode.value}) ===")
                    cmd_gain = f"{DP3_CMD} {shlex.quote(str(par_gain))}"
                    print(cmd_gain)
                    _run_to_log(cmd_gain, log_gain)

                # Step-3A.5: copy + shift instrument time window to Sun ROI
                if parmdb_shifted.exists() and overwrite.value:
                    shutil.rmtree(parmdb_shifted)

                parmdb_used = _copy_and_shift_instrument_to_target_ms(
                    inst_in=parmdb_raw,
                    target_ms=sun_roi_ms,
                    inst_out=parmdb_shifted,
                    pad_s=1.0,
                    overwrite=overwrite.value,
                )

                # Step-3B: apply shifted instrument to Sun ROI
                applycal_parset = _applycal_parset_text(sun_roi_ms, parmdb_used)
                wrote_apply = _write_text(par_apply, applycal_parset, force=overwrite.value)
                if wrote_apply:
                    print(f"[WRITE parset] {par_apply}")

                log_apply = sb_dir / "02_applycal.log"
                print(f"=== Step-3B applycal: {_sb_label(sb)} ===")
                cmd_apply = f"{DP3_CMD} {shlex.quote(str(par_apply))}"
                print(cmd_apply)
                _run_to_log(cmd_apply, log_apply)

                print(f"OK: {_sb_label(sb)} -> {sb_dir}\n")

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
    default_size=(1024, 1024),
    default_scale="1amin",
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

    Upgraded UI:
    - read SB list from plan_json_path
    - show SB + frequency in UI
    - only show SBs that already have Step-2 ROI products
    - simplified notebook usage:
          run_step4_wsclean_ui(plan_json_path=PLAN_JSON, step2_root=STEP2_ROOT, out_root=STEP4_ROOT)
    """

    import json
    import time
    import subprocess
    import shlex
    from pathlib import Path

    import ipywidgets as w
    from IPython.display import display

    # --------------------------------------------------
    # paths / payload
    # --------------------------------------------------
    plan_json_path = Path(plan_json_path)
    step2_root = Path(step2_root)
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    payload = json.load(open(plan_json_path, "r"))
    sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
    ctr_mhz = payload.get("ctr_mhz") or []

    if not sbs:
        raise ValueError("JSON must contain 'selected_sb' or 'selected_sbs'.")

    freq_map = dict(zip(sbs, ctr_mhz)) if len(ctr_mhz) == len(sbs) else {}

    # Command base (inside container)
    WSCLEAN_CMD = f"{engine} exec {bind} {shlex.quote(sif)} wsclean"

    # --------------------------------------------------
    # helpers
    # --------------------------------------------------
    def _run_to_log(cmd: str, log_path: Path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as f:
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, text=True)
            p.wait()
            if p.returncode != 0:
                raise RuntimeError(f"Command failed (see log): {log_path}\nCMD: {cmd}")

    def _roi_ms_path(sb: str) -> Path:
        tag = sb.replace(".MS", "")
        return step2_root / "ROI" / tag / f"SUN_{tag}_ROI.MS"

    def _image_prefix(sb: str) -> Path:
        tag = sb.replace(".MS", "")
        return out_root / tag / tag

    def _outputs_exist(prefix: Path) -> bool:
        pat1 = str(prefix) + "-image.fits"
        pat2 = str(prefix) + "-0000-image.fits"
        return Path(pat1).exists() or Path(pat2).exists()

    def _parse_size(text: str):
        s = text.strip().lower().replace("x", " ").replace(",", " ")
        parts = [p for p in s.split() if p]
        if len(parts) != 2:
            raise ValueError("Size must have two integers, e.g. '1024 1024'.")
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

        if weight_mode == "briggs":
            args += [f"-weight briggs {float(robust)}"]
        else:
            args += [f"-weight {weight_mode}"]

        args += [
            f"-mgain {float(mgain)}",
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

        ist = str(interval_start).strip()
        ied = str(interval_end).strip()
        if ist and ied:
            args.insert(args.index(f"-intervals-out {int(intervals_out)}"), f"-interval {int(ist)} {int(ied)}")
        elif ist or ied:
            raise ValueError("Interval start/end must be both set or both empty.")

        return f"{WSCLEAN_CMD} " + " ".join(args)

    def _sb_label(sb):
        if sb in freq_map:
            try:
                return f"{sb}   |   {float(freq_map[sb]):.1f} MHz"
            except Exception:
                return sb
        return sb

    def _available_sbs():
        good = []
        for sb in sbs:
            if _roi_ms_path(sb).exists():
                good.append(sb)
        return good

    available_sbs = _available_sbs()
    if len(available_sbs) == 0:
        raise FileNotFoundError(
            f"No available SUN_*_ROI.MS found under step2_root={step2_root}. "
            f"Please run Step-2 first, or check the step2_root path."
        )

    sb_options = [(_sb_label(sb), sb) for sb in available_sbs]

    # --------------------------------------------------
    # UI
    # --------------------------------------------------
    overwrite = w.Checkbox(value=False, description="Overwrite existing outputs")

    data_col = w.Dropdown(
        options=["CORR_NO_BEAM", "DATA"],
        value=default_data_column,
        description="Data column",
        layout=w.Layout(width="300px"),
    )

    size_box = w.Text(
        value=f"{default_size[0]} {default_size[1]}",
        description="Size (px)",
        placeholder="e.g. 1024 1024",
        layout=w.Layout(width="300px"),
    )

    scale_box = w.Text(
        value=str(default_scale),
        description="Scale",
        placeholder="e.g. 1amin / 50asec",
        layout=w.Layout(width="300px"),
    )

    weight_mode = w.Dropdown(
        options=["briggs", "natural", "uniform"],
        value=default_weight_mode,
        description="Weight",
        layout=w.Layout(width="300px"),
    )
    robust_box = w.FloatText(value=float(default_robust), description="Robust", layout=w.Layout(width="300px"))

    mgain_box = w.FloatText(value=float(default_mgain), description="mgain", layout=w.Layout(width="300px"))
    mem_box = w.IntText(value=int(default_mem_gb), description="mem (GB)", layout=w.Layout(width="300px"))
    niter_box = w.IntText(value=int(default_niter), description="niter", layout=w.Layout(width="300px"))
    amask_box = w.FloatText(value=float(default_auto_mask), description="auto-mask", layout=w.Layout(width="300px"))
    athr_box = w.FloatText(value=float(default_auto_threshold), description="auto-threshold", layout=w.Layout(width="300px"))

    intervals_out_box = w.IntText(
        value=int(default_intervals_out), description="intervals-out", layout=w.Layout(width="300px")
    )
    istart_box = w.Text(value=str(default_interval_start), description="interval start", layout=w.Layout(width="300px"))
    iend_box = w.Text(value=str(default_interval_end), description="interval end", layout=w.Layout(width="300px"))

    pol_box = w.Dropdown(options=["I"], value="I", description="Pol", layout=w.Layout(width="300px"))

    no_reorder = w.Checkbox(value=True, description="-no-reorder")
    no_update_model_required = w.Checkbox(value=True, description="-no-update-model-required")

    default_value = ("SB359.MS",) if "SB359.MS" in available_sbs else (available_sbs[0],)

    sb_select = w.SelectMultiple(
        options=sb_options,
        value=default_value,
        description="SBs",
        layout=w.Layout(width="95%", height="220px"),
    )

    run_btn = w.Button(description="Run Step-4 (WSClean)", button_style="success")
    out = w.Output()

    def _toggle_robust(*_):
        robust_box.disabled = (weight_mode.value != "briggs")

    weight_mode.observe(_toggle_robust, names="value")
    _toggle_robust()

    # --------------------------------------------------
    # run
    # --------------------------------------------------
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
            print(f"Available SBs from Step-2: {len(available_sbs)}")
            print("SBs:", ", ".join([_sb_label(sb) for sb in chosen]))
            print()

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
                    print(f"[SKIP] {_sb_label(sb)}: ROI MS not found: {ms_in}")
                    continue

                if (not overwrite.value) and _outputs_exist(prefix):
                    print(f"[SKIP] {_sb_label(sb)}: outputs exist (enable overwrite to rerun): {prefix}*")
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

                print(f"=== Step-4 WSClean: {_sb_label(sb)} ===")
                print(cmd)
                _run_to_log(cmd, log_path)
                print(f"OK: {_sb_label(sb)} -> {sb_dir}")
                print(f"Log: {log_path}")
                print()

            print(f"All done. Elapsed: {time.time() - t0:.1f} s")

    run_btn.on_click(on_run)

    # --------------------------------------------------
    # layout
    # --------------------------------------------------
    left = w.VBox([sb_select])

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
        layout=w.Layout(width="320px"),
    )

    right = w.VBox(
        [
            intervals_out_box,
            istart_box,
            iend_box,
            pol_box,
            w.HBox([no_reorder, no_update_model_required]),
        ],
        layout=w.Layout(width="340px"),
    )

    header = w.HBox([overwrite, run_btn])

    display(w.VBox([header, w.HBox([left, mid, right]), out]))



def run_step4_quicklook_ui(
    step4_root,
    out_root=None,
    plan_json_path=None,
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
      - only list SBs that actually contain image.fits
      - optionally read frequency labels from plan_json_path
      - multi-select FITS
      - generate png quicklooks (+ optional video)

    Input:
      step4_root: e.g. /data/.../step4_outputs_YYYYMMDD
      out_root:   where to save quicklooks; default step4_root/quicklook
      plan_json_path: optional selected_sb_pair_list.json for SB->freq labeling

    Output:
      out_root/SBxxx/*.png
      out_root/SBxxx/*quicklook.mp4 (optional)
      out_root/SBxxx/quicklook.log
    """
    import os
    import re
    import glob
    import json
    import shutil
    import subprocess
    import numpy as np
    import ipywidgets as w
    from pathlib import Path
    from IPython.display import display

    # --------------------------------------------------
    # helpers
    # --------------------------------------------------
    def _sb_tag(sb):
        return sb.replace(".MS", "").replace("SB", "SB")

    def _parse_tindex(fname: str):
        m = re.search(r"-t(\d+)-", os.path.basename(fname))
        return int(m.group(1)) if m else 10**9

    def _ensure_dir(p: Path):
        p.mkdir(parents=True, exist_ok=True)
        return p

    def _load_sort_fits(files):
        return sorted(files, key=lambda x: (_parse_tindex(x), os.path.basename(x)))

    def _list_fits(step4_root_path: Path, sb: str):
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
        return _load_sort_fits(sorted(set(files)))

    def _find_available_sbs(step4_root_path: Path):
        out = []
        for p in sorted(step4_root_path.glob("SB*")):
            if not (p.is_dir() and re.match(r"SB\d{3}$", p.name)):
                continue
            sb = p.name + ".MS"
            files = _list_fits(step4_root_path, sb)
            if len(files) > 0:
                out.append(sb)
        return out

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

    # --------------------------------------------------
    # plotting core (kept from your original version)
    # --------------------------------------------------
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

        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image: {fits_path} shape={getattr(data, 'shape', None)}")

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

        freq_Hz = hdr.get("CRVAL3", None)
        frequency = (freq_Hz * u.Hz) if freq_Hz is not None else (np.nan * u.Hz)

        cdelt1 = abs(hdr.get("CDELT1", np.nan)) * u.deg
        cdelt2 = abs(hdr.get("CDELT2", np.nan)) * u.deg
        cdelt1 = cdelt1.to(u.arcsec) if np.isfinite(cdelt1.value) else (np.nan * u.arcsec)
        cdelt2 = cdelt2.to(u.arcsec) if np.isfinite(cdelt2.value) else (np.nan * u.arcsec)

        rsub = None
        used_sunpy = True
        try:
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
                rotation_angle=-P1,
                wavelength=frequency.to(u.MHz) if np.isfinite(frequency.value) else None,
                observatory="NenuFAR (Nançay)",
            )

            rmap = sunpy.map.Map(data, new_header)
            rmap_rot = rmap.rotate()

            hw = float(crop_half_width_arcsec) * u.arcsec
            bl = SkyCoord(-hw, -hw, frame=rmap_rot.coordinate_frame)
            tr = SkyCoord(hw, hw, frame=rmap_rot.coordinate_frame)
            rsub = rmap_rot.submap(bl, top_right=tr)
        except Exception:
            used_sunpy = False
            rsub = None

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
                    rsub.draw_contours(np.array(contours_perc) * u.percent, colors="k",
                                       linewidths=1.1, alpha=0.8, axes=ax)
                except Exception:
                    pass
        else:
            ax = fig.add_subplot(gs[0, 0])
            cax = fig.add_subplot(gs[0, 1])

            vmin, vmax = np.nanpercentile(data, clim_pct)
            im = ax.imshow(data, origin="lower", vmin=vmin, vmax=vmax)
            ax.set_xlabel("X (pix)")
            ax.set_ylabel("Y (pix)")

        cbar = fig.colorbar(im, cax=cax)
        sf = ScalarFormatter(useMathText=True)
        sf.set_powerlimits((-2, 3))
        cbar.formatter = sf
        cbar.update_ticks()
        cbar.set_label("Stokes I (arb.)", labelpad=12)

        if draw_beam:
            try:
                from matplotlib.patches import Ellipse
                import astropy.units as u
                bmaj = hdr.get("BMAJ"); bmin = hdr.get("BMIN"); bpa = hdr.get("BPA")
                if bmaj is not None and bmin is not None and bpa is not None:
                    bmaj_as = (abs(bmaj) * u.deg).to_value(u.arcsec)
                    bmin_as = (abs(bmin) * u.deg).to_value(u.arcsec)
                    x0 = ax.get_xlim()[0] + 0.08*(ax.get_xlim()[1]-ax.get_xlim()[0])
                    y0 = ax.get_ylim()[0] + 0.08*(ax.get_ylim()[1]-ax.get_ylim()[0])
                    e = Ellipse((x0, y0), width=bmaj_as, height=bmin_as, angle=-float(bpa),
                                edgecolor="w", facecolor="none", lw=1.5)
                    ax.add_patch(e)
            except Exception:
                pass

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
        from pathlib import Path
        import os
        import shutil
        import subprocess

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
        return str(out_mp4)

    def _make_video_imageio(png_list, out_gif, fps=8):
        import imageio.v2 as imageio
        imgs = []
        for p in png_list:
            imgs.append(imageio.imread(p))
        out_gif = str(out_gif)
        imageio.mimsave(out_gif, imgs, duration=1.0/max(fps, 1))
        return out_gif

    # --------------------------------------------------
    # roots / frequency mapping
    # --------------------------------------------------
    step4_root = Path(step4_root)
    if out_root is None:
        out_root = step4_root / "quicklook"
    out_root = Path(out_root)

    freq_map = {}
    if plan_json_path is not None:
        plan_json_path = Path(plan_json_path)
        if plan_json_path.exists():
            payload = json.load(open(plan_json_path, "r"))
            sbs = payload.get("selected_sb") or payload.get("selected_sbs") or []
            ctr_mhz = payload.get("ctr_mhz") or []
            if len(sbs) == len(ctr_mhz):
                freq_map = dict(zip(sbs, ctr_mhz))

    def _sb_label(sb):
        if sb in freq_map:
            try:
                return f"{sb}   |   {float(freq_map[sb]):.1f} MHz"
            except Exception:
                return sb
        return sb

    sbs = _find_available_sbs(step4_root)
    if len(sbs) == 0:
        raise ValueError(f"No SB directories with *-image.fits found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    sb_options = [(_sb_label(sb), sb) for sb in sbs]

    # --------------------------------------------------
    # UI
    # --------------------------------------------------
    sb_dd = w.Dropdown(options=sb_options, value=default_sb_ms, description="SB", layout=w.Layout(width="340px"))
    refresh_btn = w.Button(description="Refresh FITS list", button_style="")
    fits_sel = w.SelectMultiple(options=[], description="FITS", layout=w.Layout(width="95%", height="220px"))

    overwrite_cb = w.Checkbox(value=False, description="Overwrite PNG if exists")
    make_video_cb = w.Checkbox(value=False, description="Make video (mp4/gif)")
    fps_in = w.IntText(value=8, description="FPS", layout=w.Layout(width="180px"))

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop half-width", layout=w.Layout(width="260px"))
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low %", layout=w.Layout(width="220px"))
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high %", layout=w.Layout(width="220px"))
    contours_in = w.Text(value=",".join(map(str, default_contours)), description="Contours %", layout=w.Layout(width="520px"))
    draw_beam_cb = w.Checkbox(value=True, description="Draw beam")
    run_btn = w.Button(description="Generate Quicklooks", button_style="success")

    out = w.Output()

    from astropy.coordinates import EarthLocation
    import astropy.units as u
    site = EarthLocation(lat=float(site_lat_deg)*u.deg, lon=float(site_lon_deg)*u.deg, height=float(site_height_m)*u.m)

    # --------------------------------------------------
    # callbacks
    # --------------------------------------------------
    def _refresh(_=None):
        files = _list_fits(step4_root, sb_dd.value)
        fits_sel.options = files
        fits_sel.value = tuple(files[:1]) if len(files) > 0 else ()
        with out:
            out.clear_output()
            print(f"Step4 root: {step4_root}")
            print(f"SB: {_sb_label(sb_dd.value)}  -> found {len(files)} image.fits")
            if len(files) > 0:
                print("First:", os.path.basename(files[0]))
                print("Last :", os.path.basename(files[-1]))

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

            print("Step4 root:", step4_root)
            print("SB:", _sb_label(sb))
            print("Selected:", len(files))
            print("Out:", sb_out)
            print("Overwrite:", overwrite_cb.value)
            print("Make video:", make_video_cb.value, "FPS:", fps_in.value)

            clim = (float(clim_low.value), float(clim_high.value))
            contours = _parse_contours(contours_in.value)
            crop_hw = float(crop_in.value)

            pngs = []
            for fpath in files:
                out_png = sb_out / (
                    Path(fpath).name.replace("-image.fits", "_quicklook.png").replace(".fits", "_quicklook.png")
                )
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
                mp4 = sb_out / f"{tag}_quicklook.mp4"
                gif = sb_out / f"{tag}_quicklook.gif"
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

