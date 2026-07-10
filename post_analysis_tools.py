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

from scipy.optimize import curve_fit


# ============================================================
# Shared helpers
# ============================================================

def _infer_event_tag_from_path(path_like) -> str:
    s = str(path_like)
    m = re.search(r"(20\d{6})", s)
    if m:
        return m.group(1)
    return datetime.datetime.now().strftime("%Y%m%d")


def _default_out_root(prefix: str, event_tag: str, base="/data/jzhang/nenufar_workflows") -> Path:
    return Path(base) / f"{prefix}_{event_tag}"


def _sb_tag(sb):
    return sb.replace(".MS", "").replace("SB", "SB")


def _parse_tindex(fname: str):
    m = re.search(r"-t(\d+)-", os.path.basename(fname))
    return int(m.group(1)) if m else 10**9


def _load_sort_fits(files):
    return sorted(files, key=lambda x: (_parse_tindex(x), os.path.basename(x)))


def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    return p


def _find_sbs(root: Path):
    sbs = []
    for p in sorted(root.glob("SB*")):
        if p.is_dir() and re.match(r"SB\d{3}$", p.name):
            sbs.append(p.name + ".MS")
    return sbs


def _list_step4_image_fits(step4_root: Path, sb: str):
    tag = _sb_tag(sb)
    sb_dir = step4_root / tag
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


def _list_step5b_corr_fits(step5b_root: Path, sb: str):
    if step5b_root is None:
        return []
    tag = _sb_tag(sb)
    corr_dir = step5b_root / tag / "corr_fits"
    if not corr_dir.exists():
        return []
    files = sorted([str(p) for p in corr_dir.glob("*.fits") if p.is_file()])
    return files


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


def _guess_obstime(hdr):
    for k in ["DATE-OBS", "DATEOBS", "DATE_OBS", "DATE"]:
        if k in hdr:
            try:
                return Time(hdr[k])
            except Exception:
                pass
    return Time.now()


def _build_rot_hpc_map(fpath: str, site: EarthLocation):
    data, hdr = _read_2d_data_and_header(fpath)
    obstime = _guess_obstime(hdr)

    freq_Hz = hdr.get("CRVAL3", None)
    frequency = (freq_Hz * u.Hz) if freq_Hz is not None else (np.nan * u.Hz)

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
        rotation_angle=-P1,
        wavelength=frequency.to(u.MHz) if np.isfinite(frequency.value) else None,
        observatory="NenuFAR (Nançay)",
    )

    rmap = sunpy.map.Map(data, new_header)
    rmap_rot = rmap.rotate()
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


def _solve_shift(m_sub, fit):
    if fit is None:
        raise RuntimeError("Gaussian centroid fit failed.")

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
    imageio.mimsave(out_gif, imgs, duration=1.0 / max(int(fps), 1))
    return out_gif


def run_iocorr_solve_ui(
    step4_root,
    out_root=None,
    plan_json_path=None,
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
    Ionospheric-offset correction tool A: solve quiet-Sun offset

    - scan step4_root/SBxxx/*-image.fits
    - choose ONE quiet-Sun FITS + ROI (arcsec in rotated HPC)
    - compute centroid (2D Gaussian on ROI)
    - derive offset solution relative to solar-disk center
    - save solution JSON for later apply step
    - save diagnostic preview PNG

    Output:
      out_root/SBxxx/step5a_solution.json
      out_root/SBxxx/step5a_quiet_preview.png
    """

    step4_root = Path(step4_root).expanduser()

    if out_root is None:
        event_tag = _infer_event_tag_from_path(step4_root)
        out_root = _default_out_root("step_iocorrect_outputs", event_tag)
    out_root = Path(out_root).expanduser()

    # optional frequency mapping from plan json
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

    site = EarthLocation(
        lat=float(site_lat_deg) * u.deg,
        lon=float(site_lon_deg) * u.deg,
        height=float(site_height_m) * u.m,
    )

    def _plot_quiet_preview(fpath, roi, crop_hw_as, clim_pct, out_png=None):
        m_rot, hdr, obstime, freq = _build_rot_hpc_map(fpath, site=site)
        m_sub = _submap_centered(m_rot, crop_hw_as)
        roi_map = _submap_from_roi(m_sub, roi)

        fit = _centroid_gauss(roi_map, thresh_frac=thresh_frac, min_points=min_points)

        vmin, vmax = np.nanpercentile(m_sub.data, clim_pct)

        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(1, 2, 1, projection=m_sub)
        ax2 = fig.add_subplot(1, 2, 2, projection=roi_map)

        # left panel
        m_sub.plot(axes=ax1, cmap="viridis", vmin=vmin, vmax=vmax)
        try:
            m_sub.draw_limb(axes=ax1)
            m_sub.draw_grid(axes=ax1)
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
        ax1.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0],
                 transform=ax1.get_transform("pixel"), lw=2)

        # right panel
        roi_map.plot(axes=ax2, cmap="viridis", vmin=vmin, vmax=vmax)
        try:
            roi_map.draw_limb(axes=ax2)
            roi_map.draw_grid(axes=ax2)
        except Exception:
            pass
        ax2.set_title("ROI + markers")

        # disk center
        c0 = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=roi_map.coordinate_frame)
        ax2.plot_coord(c0, marker="x", markersize=10, mew=2)

        if fit is not None:
            ax2.plot_coord(fit["cen_world"], marker="+", markersize=14, mew=2)
            ax2.text(
                0.02, 0.98,
                f'centroid = ({fit["cen_Tx_as"]:.1f}", {fit["cen_Ty_as"]:.1f}")',
                transform=ax2.transAxes,
                ha="left", va="top",
                bbox=dict(fc=(1, 1, 1, 0.6), ec="none", pad=2),
            )

        plt.tight_layout()

        if out_png is not None:
            out_png = Path(out_png)
            out_png.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(out_png, dpi=150, bbox_inches="tight")

        plt.show()
        plt.close(fig)

        return m_sub, fit, obstime, freq

    # ---------------- UI ----------------
    sbs = _find_sbs(step4_root)
    if len(sbs) == 0:
        raise ValueError(f"No SBxxx directories found under: {step4_root}")

    # only keep SBs that actually have image.fits
    sbs = [sb for sb in sbs if len(_list_step4_image_fits(step4_root, sb)) > 0]
    if len(sbs) == 0:
        raise ValueError(f"No SBs with image.fits found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    sb_dd = w.Dropdown(
        options=[(_sb_label(sb), sb) for sb in sbs],
        value=default_sb_ms,
        description="SB",
        layout=w.Layout(width="340px"),
    )
    refresh_btn = w.Button(description="Refresh FITS", button_style="")
    quiet_dd = w.Dropdown(options=[], description="Quiet FITS", layout=w.Layout(width="820px"))

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop hw (as)")
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low%")
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high%")

    xmin0, xmax0, ymin0, ymax0 = default_roi_arcsec
    xmin_w = w.FloatText(value=float(xmin0), description="xmin")
    xmax_w = w.FloatText(value=float(xmax0), description="xmax")
    ymin_w = w.FloatText(value=float(ymin0), description="ymin")
    ymax_w = w.FloatText(value=float(ymax0), description="ymax")

    run_btn = w.Button(description="Solve IO offset", button_style="warning")
    out = w.Output()

    def _refresh(_=None):
        files = _list_step4_image_fits(step4_root, sb_dd.value)
        files = _load_sort_fits(files)
        quiet_dd.options = files
        if len(files) > 0:
            quiet_dd.value = files[0]
        with out:
            clear_output()
            print("Step4 root:", step4_root)
            print("Out root  :", out_root)
            print("SB:", _sb_label(sb_dd.value), "->", len(files), "image.fits")
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

            tag = _sb_tag(sb_dd.value)
            sb_out = _ensure_dir(out_root / tag)
            preview_png = sb_out / "step5a_quiet_preview.png"

            m_sub, fit, obstime, freq = _plot_quiet_preview(
                quiet_f, roi, crop_hw, clim, out_png=preview_png
            )

            if fit is None:
                print("\n[FAIL] centroid fit failed. Try adjusting ROI or thresholds.")
                return

            dx_pix, dy_pix, cen_pix, tgt_pix = _solve_shift(m_sub, fit)

            print("\ncentroid (arcsec):", fit["cen_Tx_as"], fit["cen_Ty_as"])
            print("centroid pix (cropped map):", cen_pix)
            print("target pix (cropped map):", tgt_pix)
            print("dx_pix, dy_pix:", dx_pix, dy_pix)

            sol_path = sb_out / "step5a_solution.json"

            sol = dict(
                tool="iocorr_solve",
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
                note="dx_pix,dy_pix defined on CROPPED+ROTATED HPC map grid; "
                     "WCS correction step typically uses centroid_arcsec (Tx,Ty).",
            )

            with open(sol_path, "w") as f:
                json.dump(sol, f, indent=2)

            print("\nSaved:", sol_path)
            print("Preview:", preview_png)

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

def run_iocorr_apply_ui(
    step4_root,
    out_root=None,
    plan_json_path=None,
    default_sb="SB410",
    default_crop_half_width_arcsec=5000,
    default_clim_pct=(5, 99),
    default_contours=(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95),
    site_lat_deg=47.382,
    site_lon_deg=2.195,
    site_height_m=136.0,
):
    """
    Ionospheric-offset correction tool B: apply WCS correction to FITS

    - read step5a_solution.json (per SB)
    - apply correction by adjusting FITS WCS (CRVAL1/2 only)
    - write corrected FITS to out_root/SBxxx/corr_fits/
    - generate BEFORE/AFTER quicklooks
    - optional video mp4/gif
    """

    step4_root = Path(step4_root).expanduser()

    if out_root is None:
        event_tag = _infer_event_tag_from_path(step4_root)
        out_root = _default_out_root("step_iocorrect_outputs", event_tag)
    out_root = Path(out_root).expanduser()

    site = EarthLocation(
        lat=float(site_lat_deg) * u.deg,
        lon=float(site_lon_deg) * u.deg,
        height=float(site_height_m) * u.m,
    )

    # optional frequency mapping
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

    def _find_available_sbs():
        sbs = _find_sbs(step4_root)
        sbs = [sb for sb in sbs if len(_list_step4_image_fits(step4_root, sb)) > 0]
        return sbs

    def _load_solution_json(sb_tag: str):
        sol_path = out_root / sb_tag / "step5a_solution.json"
        if not sol_path.exists():
            raise FileNotFoundError(f"Missing step5a_solution.json: {sol_path}")
        with open(sol_path, "r") as f:
            sol = json.load(f)
        return sol, sol_path

    def _parse_centroid_arcsec(sol_dict):
        # preferred format written by run_iocorr_solve_ui
        c = sol_dict.get("centroid_arcsec", {})
        if isinstance(c, dict):
            tx = float(c.get("Tx", np.nan))
            ty = float(c.get("Ty", np.nan))
            if np.isfinite(tx) and np.isfinite(ty):
                return tx, ty

        # fallback older formats
        for kx, ky in [
            ("cen_Tx_as", "cen_Ty_as"),
            ("Tx_as", "Ty_as"),
            ("Tx", "Ty"),
        ]:
            tx = float(sol_dict.get(kx, np.nan))
            ty = float(sol_dict.get(ky, np.nan))
            if np.isfinite(tx) and np.isfinite(ty):
                return tx, ty

        raise ValueError("Cannot parse centroid arcsec from solution JSON.")

    def _apply_wcs_correction_to_header(hdr_in, Tx_as, Ty_as):
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

        ref_hpc_new = SkyCoord(
            (ref_hpc_old.Tx - Tx_as * u.arcsec),
            (ref_hpc_old.Ty - Ty_as * u.arcsec),
            frame=ref_hpc_old.frame,
        )
        ref_gcrs_new = ref_hpc_new.transform_to(ref_gcrs_old.frame)

        hdr["CRVAL1"] = ref_gcrs_new.ra.to_value(cunit1)
        hdr["CRVAL2"] = ref_gcrs_new.dec.to_value(cunit2)
        hdr.add_history(
            f"IOCORR APPLY: CRVAL shifted using centroid offset Tx={Tx_as:.6f} arcsec, Ty={Ty_as:.6f} arcsec"
        )
        return hdr

    def _build_hpc_rot_submap_from_datahdr(data2d, hdr, crop_half_width_arcsec):
        data = data2d
        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image, shape={getattr(data, 'shape', None)}")

        obstime = _guess_obstime(hdr)

        freq_Hz = hdr.get("CRVAL3", None)
        frequency = (freq_Hz * u.Hz) if freq_Hz is not None else (np.nan * u.Hz)

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
            rotation_angle=-P1,
            wavelength=frequency.to(u.MHz) if np.isfinite(frequency.value) else None,
            observatory="NenuFAR (Nançay)",
        )

        rmap = sunpy.map.Map(data, new_header)
        rmap_rot = rmap.rotate()

        hw = float(crop_half_width_arcsec) * u.arcsec
        bl = SkyCoord(-hw, -hw, frame=rmap_rot.coordinate_frame)
        tr = SkyCoord(hw, hw, frame=rmap_rot.coordinate_frame)
        try:
            rsub = rmap_rot.submap(bl, top_right=tr)
        except Exception:
            rsub = rmap_rot

        return rsub, obstime, frequency

    def _quicklook_before_after(
        fits_path,
        out_png,
        crop_half_width_arcsec,
        clim_pct,
        contours_perc,
        hdr_after=None,
        overwrite=False,
        cmap="viridis",
        draw_beam=True,
    ):
        out_png = Path(out_png)
        if (not overwrite) and out_png.exists():
            return str(out_png)

        with fits.open(fits_path) as hdul:
            hdr0 = hdul[0].header
            data0 = hdul[0].data

        data0_2d = data0
        while data0_2d is not None and getattr(data0_2d, "ndim", 0) > 2:
            data0_2d = data0_2d[0]
        data0_2d = np.squeeze(data0_2d)

        m_before, obstime, frequency = _build_hpc_rot_submap_from_datahdr(
            data0_2d, hdr0, crop_half_width_arcsec
        )

        hdr1 = hdr_after if hdr_after is not None else hdr0
        m_after, _, _ = _build_hpc_rot_submap_from_datahdr(
            data0_2d, hdr1, crop_half_width_arcsec
        )

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
                    mm.draw_contours(
                        np.array(contours_perc) * u.percent,
                        colors="k",
                        linewidths=1.1,
                        alpha=0.8,
                        axes=ax,
                    )
                except Exception:
                    pass

        axL.set_title("BEFORE (step4)", fontsize=13)
        axR.set_title("AFTER (WCS-corrected)", fontsize=13)

        cbar = fig.colorbar(imL, cax=cax)
        sf = ScalarFormatter(useMathText=True)
        sf.set_powerlimits((-2, 3))
        cbar.formatter = sf
        cbar.update_ticks()
        cbar.set_label("Stokes I (arb.)", labelpad=12)

        left_title = f"{os.path.basename(fits_path)}   {obstime.isot}"
        fig.text(0.02, 0.98, left_title, ha="left", va="top", fontsize=16)

        if np.isfinite(frequency.value):
            freq_str = f"{frequency.to_value(u.MHz):.1f} MHz"
            fig.text(
                0.98, 0.98, freq_str,
                ha="right", va="top", fontsize=15,
                bbox=dict(fc=(1, 1, 1, 0.7), ec="0.7", pad=2),
            )

        if draw_beam:
            try:
                from matplotlib.patches import Ellipse
                bmaj = hdr0.get("BMAJ")
                bmin = hdr0.get("BMIN")
                bpa = hdr0.get("BPA")
                if bmaj is not None and bmin is not None and bpa is not None:
                    bmaj_as = (abs(bmaj) * u.deg).to_value(u.arcsec)
                    bmin_as = (abs(bmin) * u.deg).to_value(u.arcsec)
                    x0 = axL.get_xlim()[0] + 0.08 * (axL.get_xlim()[1] - axL.get_xlim()[0])
                    y0 = axL.get_ylim()[0] + 0.08 * (axL.get_ylim()[1] - axL.get_ylim()[0])
                    e = Ellipse(
                        (x0, y0),
                        width=bmaj_as,
                        height=bmin_as,
                        angle=-float(bpa),
                        edgecolor="w",
                        facecolor="none",
                        lw=1.5,
                    )
                    axL.add_patch(e)
            except Exception:
                pass

        out_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return str(out_png)

    # ---------------- UI ----------------
    sbs = _find_available_sbs()
    if len(sbs) == 0:
        raise ValueError(f"No SBs with image.fits found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    sb_dd = w.Dropdown(
        options=[(_sb_label(sb), sb) for sb in sbs],
        value=default_sb_ms,
        description="SB",
        layout=w.Layout(width="340px"),
    )
    refresh_btn = w.Button(description="Refresh FITS", button_style="")
    fits_sel = w.SelectMultiple(options=[], description="FITS", layout=w.Layout(width="95%", height="220px"))

    select_all_cb = w.Checkbox(value=True, description="Select ALL image.fits")
    overwrite_cb = w.Checkbox(value=False, description="Overwrite outputs")
    make_video_cb = w.Checkbox(value=False, description="Make video (mp4/gif)")
    fps_in = w.IntText(value=8, description="FPS", layout=w.Layout(width="180px"))

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop hw (as)")
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low%")
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high%")
    contours_in = w.Text(
        value=",".join(map(str, default_contours)),
        description="Contours%",
        layout=w.Layout(width="520px"),
    )

    run_btn = w.Button(description="Apply IO correction", button_style="success")
    out = w.Output()

    def _refresh(_=None):
        files = _list_step4_image_fits(step4_root, sb_dd.value)
        files = _load_sort_fits(files)
        fits_sel.options = files
        with out:
            clear_output()
            print("Step4 root:", step4_root)
            print("Out root  :", out_root)
            print("SB:", _sb_label(sb_dd.value), "->", len(files), "image.fits")
            if len(files) > 0:
                print("First:", os.path.basename(files[0]))
                print("Last :", os.path.basename(files[-1]))

    def _on_run(_=None):
        with out:
            clear_output()

            sb = sb_dd.value
            tag = _sb_tag(sb)

            if select_all_cb.value:
                chosen = list(fits_sel.options)
            else:
                chosen = list(fits_sel.value)

            if len(chosen) == 0:
                print("No FITS selected.")
                return

            try:
                sol, sol_path = _load_solution_json(tag)
                Tx_as, Ty_as = _parse_centroid_arcsec(sol)
            except Exception as e:
                print("[FAIL] cannot load/parse solution:", e)
                return

            clim = (float(clim_low.value), float(clim_high.value))
            contours = _parse_contours(contours_in.value)
            crop_hw = float(crop_in.value)

            sb_out = _ensure_dir(out_root / tag)
            corr_dir = _ensure_dir(sb_out / "corr_fits")
            ql_dir = _ensure_dir(sb_out / "quicklook_step5b")
            log_path = sb_out / "step5b_apply.log"

            print("SB:", _sb_label(sb))
            print("Using solution:", sol_path)
            print(f"centroid offset (arcsec): Tx={Tx_as:.3f}, Ty={Ty_as:.3f}")
            print("Selected FITS:", len(chosen), "(select_all =", select_all_cb.value, ")")
            print("corr_fits  :", corr_dir)
            print("quicklooks :", ql_dir)
            print("Overwrite  :", overwrite_cb.value)
            print("Video      :", make_video_cb.value, "FPS:", int(fps_in.value))

            pngs = []
            ok = 0

            for fpath in chosen:
                fin = Path(fpath)
                fout = corr_dir / fin.name.replace("-image.fits", "-image_corrWCS.fits")
                out_png = ql_dir / fin.name.replace("-image.fits", "_step5b_before_after.png")

                try:
                    with fits.open(fin) as hdul:
                        hdr0 = hdul[0].header
                        data0 = hdul[0].data

                    hdr1 = _apply_wcs_correction_to_header(hdr0, Tx_as=Tx_as, Ty_as=Ty_as)

                    if overwrite_cb.value or (not fout.exists()):
                        fits.PrimaryHDU(data=data0, header=hdr1).writeto(fout, overwrite=True)

                    png = _quicklook_before_after(
                        fits_path=str(fin),
                        out_png=str(out_png),
                        crop_half_width_arcsec=crop_hw,
                        clim_pct=clim,
                        contours_perc=contours,
                        hdr_after=hdr1,
                        overwrite=overwrite_cb.value,
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
    sb_dd.observe(_refresh, names="value")
    run_btn.on_click(_on_run)

    _refresh()

    display(
        w.VBox([
            w.HBox([sb_dd, refresh_btn, select_all_cb, overwrite_cb, make_video_cb, fps_in]),
            fits_sel,
            w.HBox([crop_in, clim_low, clim_high]),
            w.HBox([contours_in]),
            run_btn,
            out,
        ])
    )



def run_iocorr_centroid_ui(
    step4_root,
    step5b_root=None,   # optional: corrected FITS root
    out_root=None,      # default: /data/.../step_iocentroid_outputs_YYYYMMDD
    plan_json_path=None,
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
    Ionospheric-offset correction tool C: centroid measurement

    - choose FITS from:
        * Step4 raw image.fits
        * or Step5B corrected corr_fits/*.fits
    - define ROI and fit centroid via 2D Gaussian
    - estimate centroid uncertainties following a Kontar-style approximation:
          err_x ~ sqrt(2/pi * sigma_x/sigma_y) * (deltaF/S0) * h
          err_y ~ sqrt(2/pi * sigma_y/sigma_x) * (deltaF/S0) * h
      where:
          S0     = fitted peak amplitude
          deltaF = background RMS estimated outside ROI
          h      = angular resolution ~ sqrt(BMAJ * BMIN)
    - save quicklook(s), CSV, JSONL
    - optional movie generation
    - optional beam ellipse drawing
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
    from matplotlib.patches import Ellipse

    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.io import fits
    from astropy.time import Time

    import sunpy.map
    from sunpy.coordinates import frames, sun

    from scipy.optimize import curve_fit

    # --------------------------------------------------
    # local helpers
    # --------------------------------------------------
    def _infer_event_tag_from_path(path_like) -> str:
        s = str(path_like)
        m = re.search(r"(20\d{6})", s)
        if m:
            return m.group(1)
        return datetime.datetime.now().strftime("%Y%m%d")

    def _default_out_root(prefix: str, event_tag: str) -> Path:
        return Path("/data/jzhang/nenufar_workflows") / f"{prefix}_{event_tag}"

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

    # --------------------------------------------------
    # paths
    # --------------------------------------------------
    step4_root = Path(step4_root).expanduser()
    if step5b_root is not None:
        step5b_root = Path(step5b_root).expanduser()

    if out_root is None:
        event_tag = _infer_event_tag_from_path(step4_root)
        out_root = _default_out_root("step_iocentroid_outputs", event_tag)
    out_root = Path(out_root).expanduser()

    # --------------------------------------------------
    # optional frequency mapping
    # --------------------------------------------------
    freq_map = {}
    if plan_json_path is not None:
        plan_json_path = Path(plan_json_path)
        if plan_json_path.exists():
            payload = json.load(open(plan_json_path, "r"))
            sbs0 = payload.get("selected_sb") or payload.get("selected_sbs") or []
            ctr_mhz = payload.get("ctr_mhz") or []
            if len(sbs0) == len(ctr_mhz):
                freq_map = dict(zip(sbs0, ctr_mhz))

    def _sb_label(sb):
        if sb in freq_map:
            try:
                return f"{sb}   |   {float(freq_map[sb]):.1f} MHz"
            except Exception:
                return sb
        return sb

    # --------------------------------------------------
    # site
    # --------------------------------------------------
    site = EarthLocation(
        lat=float(site_lat_deg) * u.deg,
        lon=float(site_lon_deg) * u.deg,
        height=float(site_height_m) * u.m,
    )

    # --------------------------------------------------
    # FITS -> rotated HPC map
    # --------------------------------------------------
    def _read_2d(fits_path: str):
        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image: {fits_path} shape={getattr(data, 'shape', None)}")
        return data.astype(float), hdr

    def _build_hpc_rotated_map(data, hdr):
        obstime = _guess_obstime(hdr)
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

    def _submap_from_roi_tuple(m, roi):
        xmin, xmax, ymin, ymax = roi
        p0 = SkyCoord(float(xmin) * u.arcsec, float(ymin) * u.arcsec, frame=m.coordinate_frame)
        p1 = SkyCoord(float(xmax) * u.arcsec, float(ymax) * u.arcsec, frame=m.coordinate_frame)
        return m.submap(p0, top_right=p1)

    # --------------------------------------------------
    # centroid fit
    # --------------------------------------------------
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
            popt, _pcov = curve_fit(gauss2d, (x, y), z, p0=p0, bounds=bounds, maxfev=20000)
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
            amp=float(A),
            sx_pix=float(sx_fit),
            sy_pix=float(sy_fit),
            theta_rad=float(theta_fit),
            popt=[float(v) for v in popt],
        )

    # --------------------------------------------------
    # Kontar-style error analysis
    # --------------------------------------------------
    def _estimate_background_rms_fast(m, roi):
        """
        Safer background RMS estimator using world coords sampled along
        the x and y axes separately. No unpacking of pixel_to_world().
        """
        data = np.array(m.data, dtype=float)

        ny, nx = data.shape
        xpix = np.arange(nx) * u.pix
        ypix = np.arange(ny) * u.pix

        world_x = m.pixel_to_world(xpix, np.full(nx, ny // 2) * u.pix)
        world_y = m.pixel_to_world(np.full(ny, nx // 2) * u.pix, ypix)

        try:
            x_as = world_x.Tx.to_value(u.arcsec)
        except Exception:
            x_as = world_x.spherical.lon.to_value(u.arcsec)

        try:
            y_as = world_y.Ty.to_value(u.arcsec)
        except Exception:
            y_as = world_y.spherical.lat.to_value(u.arcsec)

        xmin, xmax, ymin, ymax = roi
        x_mask = (x_as >= xmin) & (x_as <= xmax)
        y_mask = (y_as >= ymin) & (y_as <= ymax)

        roi_mask = np.outer(y_mask, x_mask)
        bg = data[~roi_mask]
        bg = bg[np.isfinite(bg)]

        if bg.size < 20:
            return np.nan

        med = np.median(bg)
        mad = np.median(np.abs(bg - med))
        return 1.4826 * mad

    def _estimate_centroid_errors_from_kontar(m, hdr, roi, fit):
        if fit is None:
            return dict(
                deltaF=np.nan,
                S0=np.nan,
                sigma_x_as=np.nan,
                sigma_y_as=np.nan,
                beam_bmaj_as=np.nan,
                beam_bmin_as=np.nan,
                h_as=np.nan,
                err_x_as=np.nan,
                err_y_as=np.nan,
            )

        S0 = float(fit["amp"])
        sx_pix = float(fit["sx_pix"])
        sy_pix = float(fit["sy_pix"])

        cdelt1_as = (abs(hdr.get("CDELT1", np.nan)) * u.deg).to_value(u.arcsec)
        cdelt2_as = (abs(hdr.get("CDELT2", np.nan)) * u.deg).to_value(u.arcsec)

        sigma_x_as = sx_pix * cdelt1_as if np.isfinite(cdelt1_as) else np.nan
        sigma_y_as = sy_pix * cdelt2_as if np.isfinite(cdelt2_as) else np.nan

        bmaj = hdr.get("BMAJ", np.nan)
        bmin = hdr.get("BMIN", np.nan)

        beam_bmaj_as = (abs(bmaj) * u.deg).to_value(u.arcsec) if np.isfinite(bmaj) else np.nan
        beam_bmin_as = (abs(bmin) * u.deg).to_value(u.arcsec) if np.isfinite(bmin) else np.nan

        if np.isfinite(beam_bmaj_as) and np.isfinite(beam_bmin_as) and beam_bmaj_as > 0 and beam_bmin_as > 0:
            h_as = np.sqrt(beam_bmaj_as * beam_bmin_as)
        else:
            h_as = np.nan

        deltaF = _estimate_background_rms_fast(m, roi)

        if (
            (not np.isfinite(S0)) or S0 <= 0 or
            (not np.isfinite(deltaF)) or deltaF < 0 or
            (not np.isfinite(sigma_x_as)) or sigma_x_as <= 0 or
            (not np.isfinite(sigma_y_as)) or sigma_y_as <= 0 or
            (not np.isfinite(h_as)) or h_as <= 0
        ):
            err_x_as = np.nan
            err_y_as = np.nan
        else:
            err_x_as = np.sqrt((2.0 / np.pi) * (sigma_x_as / sigma_y_as)) * (deltaF / S0) * h_as
            err_y_as = np.sqrt((2.0 / np.pi) * (sigma_y_as / sigma_x_as)) * (deltaF / S0) * h_as

        return dict(
            deltaF=float(deltaF) if np.isfinite(deltaF) else np.nan,
            S0=float(S0) if np.isfinite(S0) else np.nan,
            sigma_x_as=float(sigma_x_as) if np.isfinite(sigma_x_as) else np.nan,
            sigma_y_as=float(sigma_y_as) if np.isfinite(sigma_y_as) else np.nan,
            beam_bmaj_as=float(beam_bmaj_as) if np.isfinite(beam_bmaj_as) else np.nan,
            beam_bmin_as=float(beam_bmin_as) if np.isfinite(beam_bmin_as) else np.nan,
            h_as=float(h_as) if np.isfinite(h_as) else np.nan,
            err_x_as=float(err_x_as) if np.isfinite(err_x_as) else np.nan,
            err_y_as=float(err_y_as) if np.isfinite(err_y_as) else np.nan,
        )

    # --------------------------------------------------
    # draw ROI + beam
    # --------------------------------------------------
    def _draw_roi_box(ax, m, roi):
        xmin, xmax, ymin, ymax = roi
        p0 = SkyCoord(float(xmin) * u.arcsec, float(ymin) * u.arcsec, frame=m.coordinate_frame)
        p1 = SkyCoord(float(xmax) * u.arcsec, float(ymax) * u.arcsec, frame=m.coordinate_frame)
        x0, y0 = m.world_to_pixel(p0)
        x1, y1 = m.world_to_pixel(p1)
        x0 = float(np.atleast_1d(getattr(x0, "value", x0))[0])
        y0 = float(np.atleast_1d(getattr(y0, "value", y0))[0])
        x1 = float(np.atleast_1d(getattr(x1, "value", x1))[0])
        y1 = float(np.atleast_1d(getattr(y1, "value", y1))[0])
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0],
                transform=ax.get_transform("pixel"), lw=2)

    def _draw_beam(ax, m, hdr):
        try:
            bmaj = hdr.get("BMAJ", None)
            bmin = hdr.get("BMIN", None)
            bpa = hdr.get("BPA", None)
            if bmaj is None or bmin is None or bpa is None:
                return

            bmaj_as = (abs(bmaj) * u.deg).to_value(u.arcsec)
            bmin_as = (abs(bmin) * u.deg).to_value(u.arcsec)

            tx0, tx1 = -4000.0, 4000.0
            ty0, ty1 = -4000.0, 4000.0

            try:
                bl = m.bottom_left_coord
                tr = m.top_right_coord
                tx0 = float(bl.Tx.to_value(u.arcsec))
                ty0 = float(bl.Ty.to_value(u.arcsec))
                tx1 = float(tr.Tx.to_value(u.arcsec))
                ty1 = float(tr.Ty.to_value(u.arcsec))
            except Exception:
                pass

            tx_beam = tx1 - 0.12 * (tx1 - tx0)
            ty_beam = ty0 + 0.10 * (ty1 - ty0)

            beam_center = SkyCoord(tx_beam * u.arcsec, ty_beam * u.arcsec, frame=m.coordinate_frame)

            xpix, ypix = m.world_to_pixel(beam_center)
            xpix = float(np.atleast_1d(getattr(xpix, "value", xpix))[0])
            ypix = float(np.atleast_1d(getattr(ypix, "value", ypix))[0])

            cdelt1_as = (abs(hdr.get("CDELT1", np.nan)) * u.deg).to_value(u.arcsec)
            cdelt2_as = (abs(hdr.get("CDELT2", np.nan)) * u.deg).to_value(u.arcsec)

            if not np.isfinite(cdelt1_as) or not np.isfinite(cdelt2_as) or cdelt1_as == 0 or cdelt2_as == 0:
                return

            width_pix = bmaj_as / cdelt1_as
            height_pix = bmin_as / cdelt2_as

            e = Ellipse(
                (xpix, ypix),
                width=width_pix,
                height=height_pix,
                angle=-float(bpa),
                edgecolor="white",
                facecolor="none",
                lw=2.0,
                alpha=1.0,
                transform=ax.get_transform("pixel"),
                zorder=20,
            )
            ax.add_patch(e)

        except Exception as e:
            print(f"[WARN] draw_beam failed: {e}")

    # --------------------------------------------------
    # plotting
    # --------------------------------------------------
    def _quicklook_centroid_one(
        fits_path,
        out_png,
        crop_hw_as,
        clim_pct,
        contours,
        roi,
        thresh_frac,
        min_points,
        overwrite=False,
        show_beam=True,
    ):
        out_png = Path(out_png)
        if (not overwrite) and out_png.exists():
            return str(out_png), None

        data, hdr = _read_2d(fits_path)
        m0, obstime, freq_mhz = _build_hpc_rotated_map(data, hdr)
        m = _crop_center(m0, crop_hw_as)
        roi_m = _submap_from_roi_tuple(m, roi)

        fit = _centroid_gauss(roi_m, thresh_frac=thresh_frac, min_points=min_points)
        err = _estimate_centroid_errors_from_kontar(m, hdr, roi, fit)

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
                m.draw_contours(np.array(contours) * u.percent, colors="k",
                                linewidths=1.1, alpha=0.8, axes=ax)
            except Exception:
                pass

        _draw_roi_box(ax, m, roi)

        # disk center
        c0 = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=m.coordinate_frame)
        ax.plot_coord(c0, marker="x", markersize=10, mew=2)

        # centroid
        if fit is not None:
            ax.plot_coord(fit["cen_world"], marker="+", markersize=14, mew=2)

            txt = f'cen=({fit["cen_Tx_as"]:.1f}", {fit["cen_Ty_as"]:.1f}")'
            if np.isfinite(err["err_x_as"]) and np.isfinite(err["err_y_as"]):
                txt += f'\nerr=(±{err["err_x_as"]:.1f}", ±{err["err_y_as"]:.1f}")'

            ax.text(
                0.02, 0.02,
                txt,
                transform=ax.transAxes, ha="left", va="bottom",
                bbox=dict(fc=(1, 1, 1, 0.6), ec="none", pad=2),
            )

        # beam
        if show_beam:
            _draw_beam(ax, m, hdr)

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
            fig.text(
                0.98, 0.98, f"{freq_mhz.to_value(u.MHz):.1f} MHz",
                ha="right", va="top", fontsize=15,
                bbox=dict(fc=(1, 1, 1, 0.7), ec="0.7", pad=2),
            )

        out_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)

        meta = dict(
            obstime_isot=str(obstime.isot),
            freq_mhz=float(freq_mhz.to_value(u.MHz)) if np.isfinite(freq_mhz.value) else np.nan,
            cen_Tx_as=(fit["cen_Tx_as"] if fit else np.nan),
            cen_Ty_as=(fit["cen_Ty_as"] if fit else np.nan),
            amp=(fit["amp"] if fit else np.nan),
            sx_pix=(fit["sx_pix"] if fit else np.nan),
            sy_pix=(fit["sy_pix"] if fit else np.nan),
            theta_rad=(fit["theta_rad"] if fit else np.nan),
            deltaF=err["deltaF"],
            S0=err["S0"],
            sigma_x_as=err["sigma_x_as"],
            sigma_y_as=err["sigma_y_as"],
            beam_bmaj_as=err["beam_bmaj_as"],
            beam_bmin_as=err["beam_bmin_as"],
            h_as=err["h_as"],
            err_x_as=err["err_x_as"],
            err_y_as=err["err_y_as"],
            ok=(fit is not None),
        )

        return str(out_png), meta

    # --------------------------------------------------
    # movie helpers
    # --------------------------------------------------
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
        imageio.mimsave(out_gif, imgs, duration=1.0 / max(int(fps), 1))
        return out_gif

    # --------------------------------------------------
    # UI
    # --------------------------------------------------
    sbs = _find_sbs(step4_root)
    if len(sbs) == 0:
        raise ValueError(f"No SBxxx directories found under: {step4_root}")

    default_sb_ms = (default_sb + ".MS") if (default_sb and not default_sb.endswith(".MS")) else default_sb
    if default_sb_ms not in sbs:
        default_sb_ms = sbs[0]

    source_dd = w.Dropdown(
        options=["Step4 raw", "Step5B corrected"],
        value="Step4 raw",
        description="Source",
        layout=w.Layout(width="240px"),
    )
    sb_dd = w.Dropdown(
        options=[(_sb_label(sb), sb) for sb in sbs],
        value=default_sb_ms,
        description="SB",
        layout=w.Layout(width="340px"),
    )
    refresh_btn = w.Button(description="Refresh FITS", button_style="")

    fits_sel = w.SelectMultiple(
        options=[],
        description="FITS",
        layout=w.Layout(width="95%", height="220px"),
    )
    select_all_cb = w.Checkbox(value=False, description="Select ALL image.fits")
    overwrite_cb = w.Checkbox(value=False, description="Overwrite outputs")

    make_video_cb = w.Checkbox(value=False, description="Make video (mp4/gif)")
    fps_in = w.IntText(value=8, description="FPS", layout=w.Layout(width="170px"))

    crop_in = w.FloatText(value=float(default_crop_half_width_arcsec), description="Crop hw (as)")
    clim_low = w.FloatText(value=float(default_clim_pct[0]), description="CLim low%")
    clim_high = w.FloatText(value=float(default_clim_pct[1]), description="CLim high%")
    contours_in = w.Text(
        value=",".join(map(str, default_contours)),
        description="Contours%",
        layout=w.Layout(width="520px"),
    )

    xmin_in = w.FloatText(value=float(default_roi[0]), description="xmin")
    xmax_in = w.FloatText(value=float(default_roi[1]), description="xmax")
    ymin_in = w.FloatText(value=float(default_roi[2]), description="ymin")
    ymax_in = w.FloatText(value=float(default_roi[3]), description="ymax")

    thresh_in = w.FloatText(value=float(default_thresh_frac), description="thresh_frac")
    minpts_in = w.IntText(value=int(default_min_points), description="min_points")

    show_beam_cb = w.Checkbox(value=True, description="Show beam")
    show_inline_cb = w.Checkbox(value=True, description="Show images inline")
    max_inline_in = w.IntText(value=12, description="Max images", layout=w.Layout(width="170px"))

    run_btn = w.Button(description="Run centroid", button_style="success")
    out = w.Output()

    # --------------------------------------------------
    # refresh
    # --------------------------------------------------
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
            print("SB:", _sb_label(sb))
            print("Step4 root :", step4_root)
            print("Step5b root:", step5b_root if step5b_root is not None else "(None)")
            print("Found:", len(files), "fits")
            if len(files) > 0:
                print("First:", os.path.basename(files[0]))
                print("Last :", os.path.basename(files[-1]))

    # --------------------------------------------------
    # run
    # --------------------------------------------------
    def _on_run(_=None):
        with out:
            clear_output()

            sb = sb_dd.value
            tag = _sb_tag(sb)
            source = source_dd.value

            files = list(fits_sel.options) if select_all_cb.value else list(fits_sel.value)
            if len(files) == 0:
                print("No FITS selected.")
                return

            crop_hw = float(crop_in.value)
            clim = (float(clim_low.value), float(clim_high.value))
            contours = _parse_contours(contours_in.value)
            roi = (
                float(xmin_in.value),
                float(xmax_in.value),
                float(ymin_in.value),
                float(ymax_in.value),
            )
            thresh_frac = float(thresh_in.value)
            min_points = int(minpts_in.value)

            sb_out = _ensure_dir(out_root / tag)
            qdir = _ensure_dir(sb_out / "quicklook_centroid")
            log_path = sb_out / "step5c_centroid.log"
            csv_path = sb_out / "centroid_results.csv"
            jsonl_path = sb_out / "centroid_results.jsonl"

            print("Source:", source)
            print("SB:", _sb_label(sb))
            print("Out:", sb_out)
            print("ROI (as):", roi)
            print("Fit params:", f"thresh_frac={thresh_frac}", f"min_points={min_points}")
            print("Show beam:", show_beam_cb.value)
            print("----")

            if overwrite_cb.value or (not csv_path.exists()):
                with open(csv_path, "w") as f:
                    f.write(
                        "event_tag,sb,source,fname,tindex,obstime_isot,freq_mhz,"
                        "roi_xmin,roi_xmax,roi_ymin,roi_ymax,thresh_frac,min_points,"
                        "cen_Tx_as,cen_Ty_as,"
                        "amp,sx_pix,sy_pix,theta_rad,"
                        "deltaF,S0,sigma_x_as,sigma_y_as,beam_bmaj_as,beam_bmin_as,h_as,err_x_as,err_y_as,"
                        "ok,png\n"
                    )
            if overwrite_cb.value and jsonl_path.exists():
                jsonl_path.unlink()

            pngs = []
            n_ok = 0
            event_tag = _infer_event_tag_from_path(step4_root)

            for i, fpath in enumerate(files):
                fname = os.path.basename(fpath)
                tindex = _parse_tindex(fname)
                out_png = qdir / fname.replace(".fits", "_centroid.png")

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
                        show_beam=show_beam_cb.value,
                    )

                    if meta is None:
                        png, meta = _quicklook_centroid_one(
                            fits_path=fpath,
                            out_png=str(out_png),
                            crop_hw_as=crop_hw,
                            clim_pct=clim,
                            contours=contours,
                            roi=roi,
                            thresh_frac=thresh_frac,
                            min_points=min_points,
                            overwrite=True,
                            show_beam=show_beam_cb.value,
                        )

                    pngs.append(png)
                    ok = bool(meta["ok"])
                    n_ok += 1 if ok else 0

                    with open(csv_path, "a") as f:
                        f.write(
                            f"{event_tag},{tag},{source},{fname},{tindex},{meta['obstime_isot']},{meta['freq_mhz']},"
                            f"{roi[0]},{roi[1]},{roi[2]},{roi[3]},"
                            f"{thresh_frac},{min_points},"
                            f"{meta['cen_Tx_as']},{meta['cen_Ty_as']},"
                            f"{meta['amp']},{meta['sx_pix']},{meta['sy_pix']},{meta['theta_rad']},"
                            f"{meta['deltaF']},{meta['S0']},{meta['sigma_x_as']},{meta['sigma_y_as']},"
                            f"{meta['beam_bmaj_as']},{meta['beam_bmin_as']},{meta['h_as']},{meta['err_x_as']},{meta['err_y_as']},"
                            f"{int(ok)},{png}\n"
                        )

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
                        amp=meta["amp"],
                        sx_pix=meta["sx_pix"],
                        sy_pix=meta["sy_pix"],
                        theta_rad=meta["theta_rad"],
                        deltaF=meta["deltaF"],
                        S0=meta["S0"],
                        sigma_x_as=meta["sigma_x_as"],
                        sigma_y_as=meta["sigma_y_as"],
                        beam_bmaj_as=meta["beam_bmaj_as"],
                        beam_bmin_as=meta["beam_bmin_as"],
                        h_as=meta["h_as"],
                        err_x_as=meta["err_x_as"],
                        err_y_as=meta["err_y_as"],
                        ok=ok,
                        png=str(png),
                    )
                    with open(jsonl_path, "a") as f:
                        f.write(json.dumps(rec) + "\n")

                    with open(log_path, "a") as f:
                        f.write(f"OK   {fpath} -> {png}\n")

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
            w.HBox([thresh_in, minpts_in, show_beam_cb, show_inline_cb, max_inline_in]),
            run_btn,
            out,
        ])
    )

    
def run_beam_propagation_ui(
    step4_root,
    step5b_root=None,
    out_root=None,
    default_source="Step4 raw",
    default_root=None,
    default_roi=(-1500, 1000, 0, 1500),   # measurement ROI: xmin,xmax,ymin,ymax arcsec
    default_fov=(-3000, 1000, -3000, 1000),  # display FOV: xmin,xmax,ymin,ymax arcsec
    default_thresh_frac=0.5,
    default_min_points=30,
    default_projection_mode="limb",
    default_projection_angle_deg=0.0,
    default_kinematics_mode="solar altitude",
    default_background_mode="blank",
    default_aia_wavelength=193,
    default_aia_time_mode="first",
    default_aia_time_manual="",
    default_jsoc_email="",
    site_lat_deg=47.382,
    site_lon_deg=2.195,
    site_height_m=136.0,
):
    """
    Step-6 Type-III beam propagation UI (v2):
      - choose FITS across multiple SBs
      - centroid + 50% contour measurement in ROI
      - display overplot on either blank background or AIA background in independent FOV
      - kinematics mode:
          1) projected displacement
          2) solar altitude
    """
    import os
    import re
    import json
    import datetime
    from pathlib import Path

    import numpy as np
    import ipywidgets as w
    from IPython.display import display, clear_output, Image as IPyImage

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.io import fits
    from astropy.time import Time
    from astropy.constants import R_sun
    from astropy.visualization import ImageNormalize, AsinhStretch

    import sunpy.map
    from sunpy.coordinates import frames, sun
    from sunpy.net import Fido
    from sunpy.net import attrs as a

    from scipy.optimize import curve_fit

    from astropy.convolution import convolve, Gaussian2DKernel
    from reproject import reproject_interp

    C_KM_S = 299792.458
    RSUN_MM = 695.7

    # --------------------------------------------------
    # helpers: paths / tags
    # --------------------------------------------------
    def _infer_event_tag_from_any_path(p: Path) -> str:
        s = str(p)
        m = re.search(r"(20\d{6})", s)
        if m:
            return m.group(1)
        return datetime.datetime.now().strftime("%Y%m%d")

    def _default_out_root(prefix: str, event_tag: str) -> Path:
        return Path("/data/jzhang/nenufar_workflows") / f"{prefix}_{event_tag}"

    def _ensure_dir(p: Path):
        p.mkdir(parents=True, exist_ok=True)
        return p

    def _parse_tindex(fname: str):
        m = re.search(r"-t(\d+)-", os.path.basename(fname))
        return int(m.group(1)) if m else 10**9

    def _infer_sb_from_path(fpath: str):
        m = re.search(r"(SB\d{3})", str(fpath))
        return m.group(1) if m else ""

    def _load_sort_fits(files):
        return sorted(files, key=lambda x: (_parse_tindex(x), os.path.basename(x)))

    def _scan_all_image_fits(root_path: Path):
        root_path = Path(root_path)
        if not root_path.exists():
            return []

        files = []
        for p in root_path.rglob("*.fits"):
            if not p.is_file():
                continue

            base = p.name.lower()
            is_science_image = (
                base.endswith("-image.fits")
                or base.endswith("_image.fits")
                or "image_corrwcs.fits" in base
                or base.endswith("-image_corrwcs.fits")
                or base.endswith("_image_corrwcs.fits")
            )
            if not is_science_image:
                continue

            bad_keys = ["dirty", "psf", "model", "residual"]
            if any(k in base for k in bad_keys):
                continue

            files.append(str(p))

        return _load_sort_fits(sorted(set(files)))

    # --------------------------------------------------
    # helpers: FITS -> rotated HPC map
    # --------------------------------------------------
    site = EarthLocation(
        lat=float(site_lat_deg) * u.deg,
        lon=float(site_lon_deg) * u.deg,
        height=float(site_height_m) * u.m,
    )

    def _read_2d(fits_path: str):
        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
        while data is not None and getattr(data, "ndim", 0) > 2:
            data = data[0]
        data = np.squeeze(data)
        if data is None or data.ndim != 2:
            raise ValueError(f"Not a 2D image: {fits_path} shape={getattr(data, 'shape', None)}")
        return data.astype(float), hdr

    def _get_obstime(hdr):
        for k in ["DATE-OBS", "DATEOBS", "DATE_OBS", "DATE", "DATE_OBS"]:
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
        freq_mhz = frequency.to(u.MHz).value if np.isfinite(frequency.value) else np.nan
        return m, obstime, freq_mhz

    def _submap_from_box(m, box):
        xmin, xmax, ymin, ymax = box
        p0 = SkyCoord(float(xmin) * u.arcsec, float(ymin) * u.arcsec, frame=m.coordinate_frame)
        p1 = SkyCoord(float(xmax) * u.arcsec, float(ymax) * u.arcsec, frame=m.coordinate_frame)
        return m.submap(p0, top_right=p1)

    def _get_rsun_obs_as(m, hdr):
        try:
            rsun = m.rsun_obs.to_value(u.arcsec)
            if np.isfinite(rsun):
                return float(rsun)
        except Exception:
            pass

        for k in ["RSUN_OBS", "SOLAR_R", "RSUN"]:
            if k in hdr:
                try:
                    val = float(hdr[k])
                    if np.isfinite(val):
                        return val
                except Exception:
                    pass

        try:
            obstime = _get_obstime(hdr)
            ang = np.arcsin((R_sun / sun.earth_distance(obstime)).decompose().value) * u.rad
            return float(ang.to_value(u.arcsec))
        except Exception:
            return np.nan

    def _arcsec_to_km_factor(obstime):
        try:
            dist = sun.earth_distance(obstime).to(u.km)
            km_per_as = (dist * np.deg2rad(1.0 / 3600.0)).to_value(u.km)
            return float(km_per_as)
        except Exception:
            return 725.0

    # --------------------------------------------------
    # AIA helpers
    # --------------------------------------------------
    def _get_aia_reference_time(results, mode="first", manual_time_str=""):
        good = [r for r in results if r.get("ok", False)]
        if len(good) == 0:
            return None

        good = sorted(good, key=lambda r: r["obstime"].unix)

        if mode == "manual":
            try:
                return Time(manual_time_str)
            except Exception:
                return good[0]["obstime"]
        elif mode == "middle":
            return good[len(good) // 2]["obstime"]
        else:
            return good[0]["obstime"]

    def _find_cached_aia_map(ref_time, wavelength=193, search_dirs=None, max_dt=30*u.min):
        """
        Search existing downloaded AIA FITS files and return the nearest one in time.
        """
        if ref_time is None:
            return None, None

        if search_dirs is None:
            search_dirs = []

        cand_files = []
        for d in search_dirs:
            if d is None:
                continue
            d = Path(d)
            if not d.exists():
                continue
            cand_files.extend(sorted(d.rglob(f"*{wavelength}.image.fits")))
            cand_files.extend(sorted(d.rglob(f"*{wavelength}*fits")))

        cand_files = list(dict.fromkeys([str(p) for p in cand_files]))

        best_file = None
        best_dt = None

        for f in cand_files:
            try:
                with fits.open(f) as hdul:
                    hdr = hdul[0].header
                    t = _get_obstime(hdr)

                    ok_wave = True
                    if "WAVELNTH" in hdr:
                        try:
                            ok_wave = (int(hdr["WAVELNTH"]) == int(wavelength))
                        except Exception:
                            ok_wave = True

                    if not ok_wave:
                        continue

                    dt = abs((t - ref_time).to_value(u.s))
                    if best_dt is None or dt < best_dt:
                        best_dt = dt
                        best_file = f
            except Exception:
                continue

        if best_file is None:
            return None, None

        if best_dt is not None and best_dt > max_dt.to_value(u.s):
            return None, None

        try:
            best_map = sunpy.map.Map(best_file)
            return best_map, best_file
        except Exception:
            return None, None

    def _download_aia_cutout_map(ref_time, fov, wavelength=193, jsoc_email=None,
                                 sample_window=0.5*u.min, out_dir=None):
        """
        Download a single nearest-in-time AIA cutout around ref_time for the requested FOV.
        Requires a JSOC-registered email address.
        """
        if ref_time is None:
            raise ValueError("AIA reference time is None.")

        if jsoc_email is None or str(jsoc_email).strip() == "":
            jsoc_email = os.environ.get("JSOC_EMAIL", "").strip()

        if jsoc_email == "":
            raise RuntimeError(
                "JSOC email is required for AIA cutout download. "
                "Set JSOC_EMAIL in the environment or provide it in the UI."
            )

        xmin, xmax, ymin, ymax = fov
        bottom_left = SkyCoord(
            xmin * u.arcsec, ymin * u.arcsec,
            obstime=ref_time, observer="earth", frame="helioprojective"
        )
        top_right = SkyCoord(
            xmax * u.arcsec, ymax * u.arcsec,
            obstime=ref_time, observer="earth", frame="helioprojective"
        )

        cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)

        query = Fido.search(
            a.Time(ref_time - sample_window/2, ref_time + sample_window/2),
            a.Wavelength(wavelength * u.angstrom),
            a.jsoc.Series.aia_lev1_euv_12s,
            a.jsoc.Notify(jsoc_email),
            a.jsoc.Segment.image,
            cutout,
        )

        if len(query) == 0:
            raise RuntimeError("No AIA cutout query results returned.")

        qr = query[0]
        if "Start Time" in qr.colnames:
            times = Time(qr["Start Time"])
        elif "TIME" in qr.colnames:
            times = Time(qr["TIME"])
        else:
            nearest_idx = 0
            qr1 = qr[nearest_idx:nearest_idx+1]

        if "Start Time" in qr.colnames or "TIME" in qr.colnames:
            dt = np.abs((times - ref_time).to_value(u.s))
            nearest_idx = int(np.argmin(dt))
            qr1 = qr[nearest_idx:nearest_idx+1]

        if out_dir is None:
            out_dir = os.getcwd()

        files = Fido.fetch(qr1, path=str(Path(out_dir) / "{file}"))
        if len(files) == 0:
            raise RuntimeError("Fido.fetch returned no AIA files.")

        aia_map = sunpy.map.Map(files[0])
        return aia_map, files[0]

    # --------------------------------------------------
    # centroid + contour helpers
    # --------------------------------------------------
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
            xp = ct * (x_ - x0) + st * (y_ - y0)
            yp = -st * (x_ - x0) + ct * (y_ - y0)
            return (A * np.exp(-0.5 * ((xp / sx) ** 2 + (yp / sy) ** 2))).ravel()

        y0_idx, x0_idx = np.unravel_index(np.nanargmax(z2), z2.shape)
        s_guess = max(2.0, 0.15 * min(z2.shape))
        p0 = [amp0, float(x0_idx), float(y0_idx), 0.0, s_guess, s_guess]
        bounds = (
            [0.0, 0.0, 0.0, -np.pi, 1.0, 1.0],
            [3 * amp0, z2.shape[1] - 1, z2.shape[0] - 1, np.pi, max(z2.shape), max(z2.shape)],
        )

        try:
            popt, _ = curve_fit(gauss2d, (x, y), z, p0=p0, bounds=bounds, maxfev=20000)
        except Exception:
            return None

        A, x0_fit, y0_fit, theta_fit, sx_fit, sy_fit = popt
        cen_world = roi_map.pixel_to_world(x0_fit * u.pixel, y0_fit * u.pixel)

        fwhm_major_pix = 2.3548 * max(sx_fit, sy_fit)
        fwhm_minor_pix = 2.3548 * min(sx_fit, sy_fit)

        scale_x = abs(roi_map.scale.axis1.to_value(u.arcsec / u.pix))
        scale_y = abs(roi_map.scale.axis2.to_value(u.arcsec / u.pix))
        fwhm_major_as = fwhm_major_pix * scale_x
        fwhm_minor_as = fwhm_minor_pix * scale_y

        return dict(
            A=float(A),
            x0_pix=float(x0_fit),
            y0_pix=float(y0_fit),
            theta_rad=float(theta_fit),
            sx_pix=float(sx_fit),
            sy_pix=float(sy_fit),
            cen_Tx_as=float(cen_world.Tx.to_value(u.arcsec)),
            cen_Ty_as=float(cen_world.Ty.to_value(u.arcsec)),
            fwhm_major_as=float(fwhm_major_as),
            fwhm_minor_as=float(fwhm_minor_as),
        )

    def _extract_50_contour_world(roi_map):
        z2 = np.array(roi_map.data, dtype=float)
        if not np.isfinite(z2).any():
            return []
        vmax = np.nanmax(z2)
        if not np.isfinite(vmax) or vmax <= 0:
            return []

        level = 0.5 * vmax
        fig, ax = plt.subplots()
        segs = []
        try:
            cs = ax.contour(z2, levels=[level])
            for seg in cs.allsegs[0]:
                if len(seg) < 3:
                    continue
                xs = seg[:, 0]
                ys = seg[:, 1]
                world = roi_map.pixel_to_world(xs * u.pixel, ys * u.pixel)
                segs.append(
                    np.column_stack([
                        world.Tx.to_value(u.arcsec),
                        world.Ty.to_value(u.arcsec),
                    ])
                )
        finally:
            plt.close(fig)
        return segs

    def _measure_one_fits(fits_path, roi, thresh_frac, min_points):
        data, hdr = _read_2d(fits_path)
        m, obstime, freq_mhz = _build_hpc_rotated_map(data, hdr)
        roi_map = _submap_from_box(m, roi)

        fit = _centroid_gauss(roi_map, thresh_frac=thresh_frac, min_points=min_points)
        contours50 = _extract_50_contour_world(roi_map)

        rsun_obs_as = _get_rsun_obs_as(m, hdr)
        km_per_as = _arcsec_to_km_factor(obstime)

        rec = dict(
            fits=str(fits_path),
            fname=os.path.basename(fits_path),
            sb_tag=_infer_sb_from_path(fits_path),
            instrument="NenuFAR",
            obstime=obstime,
            obstime_isot=obstime.isot,
            freq_mhz=float(freq_mhz),
            rsun_obs_as=float(rsun_obs_as) if np.isfinite(rsun_obs_as) else np.nan,
            km_per_as=float(km_per_as),
            ok=False,
            contours50=contours50,
        )

        if fit is None:
            return rec

        r_as = float(np.hypot(fit["cen_Tx_as"], fit["cen_Ty_as"]))
        altitude_as = float(r_as - rsun_obs_as) if np.isfinite(rsun_obs_as) else np.nan
        altitude_km = float(altitude_as * km_per_as) if np.isfinite(altitude_as) else np.nan
        altitude_Mm = float(altitude_km / 1e3) if np.isfinite(altitude_km) else np.nan

        rec["ok"] = True
        rec.update(fit)
        rec["r_as"] = r_as
        rec["altitude_as"] = altitude_as
        rec["altitude_km"] = altitude_km
        rec["altitude_Mm"] = altitude_Mm
        return rec

    # --------------------------------------------------
    # kinematics helpers
    # --------------------------------------------------
    def _get_projection_axis(results, mode="limb", angle_deg=0.0):
        good = [r for r in results if r.get("ok", False)]
        if len(good) == 0:
            return np.array([1.0, 0.0]), np.array([0.0, 0.0])

        ref = np.array([good[0]["cen_Tx_as"], good[0]["cen_Ty_as"]], dtype=float)

        if mode == "custom":
            ang = np.deg2rad(float(angle_deg))
            axis = np.array([np.cos(ang), np.sin(ang)], dtype=float)
        else:
            axis = ref.copy()
            nrm = np.linalg.norm(axis)
            axis = np.array([1.0, 0.0], dtype=float) if nrm == 0 else axis / nrm

        return axis, ref

    def _fit_speed_from_series(tt, yy):
        tt = np.array(tt, dtype=float)
        yy = np.array(yy, dtype=float)

        good = np.isfinite(tt) & np.isfinite(yy)
        tt = tt[good]
        yy = yy[good]

        if len(tt) < 2 or np.unique(tt).size < 2:
            return dict(
                fit_speed_km_s=np.nan,
                fit_speed_c=np.nan,
                fit_speed_err_km_s=np.nan,
                fit_speed_err_c=np.nan,
                intercept=np.nan,
                pairwise_mean_km_s=np.nan,
                pairwise_std_km_s=np.nan,
                pairwise_mean_c=np.nan,
                pairwise_std_c=np.nan,
                n_good=len(tt),
                n_unique_times=int(np.unique(tt).size) if len(tt) > 0 else 0,
            )

        vv = []
        for i in range(1, len(tt)):
            dt = tt[i] - tt[i - 1]
            dy = yy[i] - yy[i - 1]
            vv.append(dy / dt if dt != 0 else np.nan)
        vv = np.array(vv, dtype=float)
        vv = vv[np.isfinite(vv)]

        pairwise_mean = float(np.nanmean(vv)) if vv.size >= 1 else np.nan
        pairwise_std = float(np.nanstd(vv, ddof=1)) if vv.size >= 2 else np.nan

        fit_speed = np.nan
        intercept = np.nan
        fit_err = np.nan

        try:
            coef = np.polyfit(tt, yy, 1)
            fit_speed = float(coef[0])
            intercept = float(coef[1])
        except Exception:
            pass

        if len(tt) >= 3:
            try:
                coef_cov, cov = np.polyfit(tt, yy, 1, cov=True)
                fit_speed = float(coef_cov[0])
                intercept = float(coef_cov[1])
                fit_err = float(np.sqrt(cov[0, 0])) if np.isfinite(cov[0, 0]) else np.nan
            except Exception:
                pass

        return dict(
            fit_speed_km_s=fit_speed,
            fit_speed_c=(fit_speed / C_KM_S if np.isfinite(fit_speed) else np.nan),
            fit_speed_err_km_s=fit_err,
            fit_speed_err_c=(fit_err / C_KM_S if np.isfinite(fit_err) else np.nan),
            intercept=intercept,
            pairwise_mean_km_s=pairwise_mean,
            pairwise_std_km_s=pairwise_std,
            pairwise_mean_c=(pairwise_mean / C_KM_S if np.isfinite(pairwise_mean) else np.nan),
            pairwise_std_c=(pairwise_std / C_KM_S if np.isfinite(pairwise_std) else np.nan),
            n_good=len(tt),
            n_unique_times=int(np.unique(tt).size),
        )

    def _add_kinematics(results, kinematics_mode="solar altitude", projection_mode="limb", angle_deg=0.0):
        good = [r for r in results if r.get("ok", False)]
        if len(good) == 0:
            return results, None

        good = sorted(good, key=lambda r: r["obstime"].unix)
        t0 = good[0]["obstime"].unix
        for r in good:
            r["dt_s"] = float(r["obstime"].unix - t0)

        if kinematics_mode == "projected displacement":
            axis, ref = _get_projection_axis(good, mode=projection_mode, angle_deg=angle_deg)

            for r in good:
                p = np.array([r["cen_Tx_as"], r["cen_Ty_as"]], dtype=float)
                s_as = np.dot(p - ref, axis)
                s_km = s_as * r.get("km_per_as", 725.0)
                r["proj_coord_as"] = float(s_as)
                r["proj_coord_km"] = float(s_km)
                r["proj_coord_Mm"] = float(s_km / 1e3)

            fit = _fit_speed_from_series(
                [r["dt_s"] for r in good],
                [r["proj_coord_km"] for r in good],
            )

            prev = None
            for r in good:
                if prev is None:
                    r["pairwise_speed_km_s"] = np.nan
                    r["pairwise_speed_c"] = np.nan
                else:
                    dt = r["dt_s"] - prev["dt_s"]
                    dy = r["proj_coord_km"] - prev["proj_coord_km"]
                    v = dy / dt if dt != 0 else np.nan
                    r["pairwise_speed_km_s"] = float(v) if np.isfinite(v) else np.nan
                    r["pairwise_speed_c"] = float(v / C_KM_S) if np.isfinite(v) else np.nan
                r["fit_speed_km_s"] = fit["fit_speed_km_s"]
                r["fit_speed_c"] = fit["fit_speed_c"]
                prev = r

            fit_meta = dict(
                kinematics_mode=kinematics_mode,
                projection_mode=projection_mode,
                axis=axis.tolist(),
                ref=ref.tolist(),
                y_label="Projected distance (Mm)",
                y_key="proj_coord_Mm",
                y_key_as="proj_coord_as",
                y_key_km="proj_coord_km",
                y_key_Mm="proj_coord_Mm",
                **fit,
            )
            return results, fit_meta

        else:
            fit = _fit_speed_from_series(
                [r["dt_s"] for r in good],
                [r["altitude_km"] for r in good],
            )

            prev = None
            for r in good:
                if prev is None:
                    r["pairwise_altitude_speed_km_s"] = np.nan
                    r["pairwise_altitude_speed_c"] = np.nan
                else:
                    dt = r["dt_s"] - prev["dt_s"]
                    dy = r["altitude_km"] - prev["altitude_km"]
                    v = dy / dt if dt != 0 else np.nan
                    r["pairwise_altitude_speed_km_s"] = float(v) if np.isfinite(v) else np.nan
                    r["pairwise_altitude_speed_c"] = float(v / C_KM_S) if np.isfinite(v) else np.nan
                r["altitude_fit_speed_km_s"] = fit["fit_speed_km_s"]
                r["altitude_fit_speed_c"] = fit["fit_speed_c"]
                prev = r

            fit_meta = dict(
                kinematics_mode=kinematics_mode,
                projection_mode=None,
                axis=None,
                ref=None,
                y_label="Solar altitude (Mm)",
                y_key="altitude_Mm",
                y_key_as="altitude_as",
                y_key_km="altitude_km",
                y_key_Mm="altitude_Mm",
                **fit,
            )
            return results, fit_meta

    # --------------------------------------------------
    # plotting helpers
    # --------------------------------------------------
    def _plot_overplot(results, fov, out_png, annotate_by="time",
                    background_mode="blank", aia_map=None):
        from matplotlib.lines import Line2D

        good = [r for r in results if r.get("ok", False)]
        if len(good) == 0:
            return

        xmin, xmax, ymin, ymax = fov

        # --------------------------------------------------
        # AIA mode: WCS axes, but NO full-map reprojection
        # --------------------------------------------------
        if background_mode == "AIA" and aia_map is not None:
            try:
                bg_map = _submap_from_box(aia_map, fov)

                fig = plt.figure(figsize=(9.0, 8.0))
                ax = fig.add_subplot(projection=bg_map)

                bg_map.plot(axes=ax)
                try:
                    bg_map.draw_grid(axes=ax, alpha=0.2)
                except Exception:
                    pass

                ax.patch.set_facecolor("black")
                fig.patch.set_facecolor("white")

                colors = plt.cm.cool(np.linspace(0.05, 0.95, max(10, len(good))))
                legend_entries = []

                for i, r in enumerate(good):
                    c = colors[i % len(colors)]

                    # draw 50% contours directly from world coords
                    for seg in r.get("contours50", []):
                        if seg is None or len(seg) < 2:
                            continue
                        seg_coord = SkyCoord(
                            seg[:, 0] * u.arcsec,
                            seg[:, 1] * u.arcsec,
                            frame=bg_map.coordinate_frame,
                        )
                        ax.plot_coord(seg_coord, color=c, lw=2.0, alpha=0.85)

                    # centroid directly in world coords
                    cen = SkyCoord(
                        r["cen_Tx_as"] * u.arcsec,
                        r["cen_Ty_as"] * u.arcsec,
                        frame=bg_map.coordinate_frame,
                    )
                    ax.plot_coord(
                        cen,
                        marker="o",
                        ms=7.5,
                        color=c,
                        mec="white",
                        mew=0.8,
                        linestyle="None",
                    )

                    inst_str = r.get("instrument", "Unknown")
                    time_str = r["obstime"].isot.split("T")[-1][:8]
                    freq_str = f'{r["freq_mhz"]:.1f} MHz' if np.isfinite(r["freq_mhz"]) else r["sb_tag"]
                    legend_entries.append((c, f"{inst_str} | {freq_str} | {time_str}"))

                proxies = [
                    Line2D([0], [0], color=c, marker='o', lw=2, markersize=8,
                        markeredgecolor='white', markeredgewidth=0.8)
                    for c, _ in legend_entries
                ]
                labels = [lbl for _, lbl in legend_entries]

                leg = ax.legend(
                    proxies,
                    labels,
                    loc="center left",
                    bbox_to_anchor=(1.02, 0.5),
                    frameon=True,
                    title="Selected burst images",
                    fontsize=9,
                    title_fontsize=10,
                    borderaxespad=0.0,
                )
                leg.get_frame().set_facecolor("white")
                leg.get_frame().set_edgecolor("black")

                ax.set_title("Step6: centroids + 50% contours (on AIA)")
                ax.set_xlabel("Tx (arcsec)")
                ax.set_ylabel("Ty (arcsec)")

                ax.set_xlim(xmin * u.arcsec, xmax * u.arcsec)
                ax.set_ylim(ymin * u.arcsec, ymax * u.arcsec)

                fig.tight_layout()
                fig.savefig(out_png, dpi=180, bbox_inches="tight", facecolor="white")
                plt.close(fig)
                return

            except Exception as e:
                print("[WARN] AIA/WCS plotting failed, fallback to blank:", e)

        # --------------------------------------------------
        # blank mode
        # --------------------------------------------------
        fig, ax = plt.subplots(figsize=(9.0, 8.0))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("#dbeeff")

        rsun_obs_as = good[0].get("rsun_obs_as", np.nan)
        if np.isfinite(rsun_obs_as):
            disk = Circle(
                (0, 0),
                rsun_obs_as,
                edgecolor="white",
                facecolor="none",
                lw=1.5,
                alpha=0.95,
                zorder=2,
            )
            ax.add_patch(disk)

        colors = plt.cm.cool(np.linspace(0.05, 0.95, max(10, len(good))))
        legend_handles = []
        legend_labels = []

        for i, r in enumerate(good):
            c = colors[i % len(colors)]

            for seg in r.get("contours50", []):
                ax.plot(seg[:, 0], seg[:, 1], color=c, lw=2.0, alpha=0.95, zorder=3)

            h, = ax.plot(
                r["cen_Tx_as"],
                r["cen_Ty_as"],
                marker="o",
                ms=7.5,
                color=c,
                mec="white",
                mew=0.8,
                linestyle="None",
                zorder=4,
            )

            inst_str = r.get("instrument", "Unknown")
            time_str = r["obstime"].isot.split("T")[-1][:8]
            freq_str = f'{r["freq_mhz"]:.1f} MHz' if np.isfinite(r["freq_mhz"]) else r["sb_tag"]
            label = f"{inst_str} | {freq_str} | {time_str}"

            legend_handles.append(h)
            legend_labels.append(label)

        ax.set_title("Step6: centroids + 50% contours (blank background)")
        ax.set_xlabel("Tx (arcsec)")
        ax.set_ylabel("Ty (arcsec)")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.15, color="gray")

        leg = ax.legend(
            legend_handles,
            legend_labels,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=True,
            title="Selected burst images",
            fontsize=9,
            title_fontsize=10,
            borderaxespad=0.0,
        )
        leg.get_frame().set_facecolor("white")
        leg.get_frame().set_edgecolor("black")

        fig.tight_layout()
        fig.savefig(out_png, dpi=180, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    def _plot_kinematics(results, fit_meta, out_png):
        good = [r for r in results if r.get("ok", False)]
        if len(good) == 0 or fit_meta is None:
            return

        good = sorted(good, key=lambda r: r.get("dt_s", np.nan))
        tt = np.array([r.get("dt_s", np.nan) for r in good], dtype=float)

        fig, ax = plt.subplots(figsize=(7.8, 5.2))

        if fit_meta["kinematics_mode"] == "solar altitude":
            yy = np.array([r.get("altitude_Mm", np.nan) for r in good], dtype=float)
            ax.plot(tt, yy, "o", ms=8, label="Measured positions")

            if np.isfinite(fit_meta.get("fit_speed_km_s", np.nan)) and np.unique(tt[np.isfinite(tt)]).size >= 2:
                tfit = np.linspace(np.nanmin(tt), np.nanmax(tt), 200)
                yfit_Mm = (fit_meta["fit_speed_km_s"] * tfit + fit_meta["intercept"]) / 1e3
                lbl = f'Linear fit: {fit_meta["fit_speed_km_s"]:.1f} km/s'
                if np.isfinite(fit_meta.get("fit_speed_err_km_s", np.nan)):
                    lbl += f' ± {fit_meta["fit_speed_err_km_s"]:.1f}'
                lbl += f' ({fit_meta["fit_speed_c"]:.4f} c)'
                ax.plot(tfit, yfit_Mm, "-", lw=2.2, label=lbl)

            ax.set_ylabel("Solar altitude (Mm)")
            ax.set_title("Step6: solar altitude vs time")

            def mm_to_rsun(y_mm):
                return 1.0 + y_mm / RSUN_MM

            def rsun_to_mm(y_rsun):
                return (y_rsun - 1.0) * RSUN_MM

            secax = ax.secondary_yaxis("right", functions=(mm_to_rsun, rsun_to_mm))
            secax.set_ylabel(r"Heliocentric distance ($R_\odot$)")

        else:
            yy = np.array([r.get("proj_coord_Mm", np.nan) for r in good], dtype=float)
            ax.plot(tt, yy, "o", ms=8, label="Measured positions")

            if np.isfinite(fit_meta.get("fit_speed_km_s", np.nan)) and np.unique(tt[np.isfinite(tt)]).size >= 2:
                tfit = np.linspace(np.nanmin(tt), np.nanmax(tt), 200)
                yfit_Mm = (fit_meta["fit_speed_km_s"] * tfit + fit_meta["intercept"]) / 1e3
                lbl = f'Linear fit: {fit_meta["fit_speed_km_s"]:.1f} km/s'
                if np.isfinite(fit_meta.get("fit_speed_err_km_s", np.nan)):
                    lbl += f' ± {fit_meta["fit_speed_err_km_s"]:.1f}'
                lbl += f' ({fit_meta["fit_speed_c"]:.4f} c)'
                ax.plot(tfit, yfit_Mm, "-", lw=2.2, label=lbl)

            ax.set_ylabel("Projected distance (Mm)")
            ax.set_title("Step6: projected displacement vs time")

        ax.set_xlabel("Time since first frame (s)")
        ax.legend()
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(out_png, dpi=180, bbox_inches="tight")
        plt.close(fig)

    def _jsonify_obj(obj):
        if isinstance(obj, dict):
            return {k: _jsonify_obj(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_jsonify_obj(v) for v in obj]
        elif isinstance(obj, tuple):
            return [_jsonify_obj(v) for v in obj]
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.floating, np.integer, np.bool_)):
            return obj.item()
        else:
            return obj

    # --------------------------------------------------
    # initial paths
    # --------------------------------------------------
    step4_root = Path(step4_root)
    step5b_root = Path(step5b_root) if step5b_root is not None else None

    if default_root is None:
        if default_source == "Step5B corrected" and step5b_root is not None:
            default_root = str(step5b_root)
        else:
            default_root = str(step4_root)

    event_tag = _infer_event_tag_from_any_path(Path(default_root))
    if out_root is None:
        out_root = _default_out_root("step6_outputs", event_tag)
    out_root = Path(out_root)
    _ensure_dir(out_root)

    # --------------------------------------------------
    # widgets
    # --------------------------------------------------
    source_dd = w.Dropdown(
        options=["Step4 raw", "Step5B corrected"],
        value=default_source,
        description="Source",
    )

    root_in = w.Text(
        value=str(default_root),
        description="Root",
        layout=w.Layout(width="900px"),
    )

    refresh_btn = w.Button(description="Refresh")
    add_btn = w.Button(description="Add ->")
    remove_btn = w.Button(description="<- Remove")
    add_all_btn = w.Button(description="Add all")
    clear_btn = w.Button(description="Clear")
    overwrite_cb = w.Checkbox(value=False, description="Overwrite outputs")

    avail_sel = w.SelectMultiple(
        options=[],
        value=(),
        description="Available FITS",
        layout=w.Layout(width="48%", height="260px"),
    )

    selected_sel = w.SelectMultiple(
        options=[],
        value=(),
        description="Selected FITS",
        layout=w.Layout(width="48%", height="260px"),
    )

    user_tag_in = w.Text(value="t3beam_run1", description="Run tag")

    # measurement ROI
    xmin_in = w.FloatText(value=float(default_roi[0]), description="ROI xmin")
    xmax_in = w.FloatText(value=float(default_roi[1]), description="ROI xmax")
    ymin_in = w.FloatText(value=float(default_roi[2]), description="ROI ymin")
    ymax_in = w.FloatText(value=float(default_roi[3]), description="ROI ymax")

    # display FOV
    fov_xmin_in = w.FloatText(value=float(default_fov[0]), description="FOV xmin")
    fov_xmax_in = w.FloatText(value=float(default_fov[1]), description="FOV xmax")
    fov_ymin_in = w.FloatText(value=float(default_fov[2]), description="FOV ymin")
    fov_ymax_in = w.FloatText(value=float(default_fov[3]), description="FOV ymax")

    thresh_in = w.FloatText(value=float(default_thresh_frac), description="thresh_frac")
    minpts_in = w.IntText(value=int(default_min_points), description="min_points")

    bg_mode_dd = w.Dropdown(
        options=["blank", "AIA"],
        value=default_background_mode,
        description="Background",
    )
    aia_wave_in = w.IntText(value=int(default_aia_wavelength), description="AIA Å")
    aia_time_mode_dd = w.Dropdown(
        options=["first", "middle", "manual"],
        value=default_aia_time_mode,
        description="AIA time",
    )
    aia_time_manual_in = w.Text(
        value=str(default_aia_time_manual),
        description="Manual time",
        placeholder="YYYY-MM-DDTHH:MM:SS",
    )
    jsoc_email_in = w.Text(
        value=str(default_jsoc_email) if str(default_jsoc_email).strip() != "" else os.environ.get("JSOC_EMAIL", ""),
        description="JSOC email",
        placeholder="registered_email@example.com",
    )

    kin_mode_dd = w.Dropdown(
        options=["projected displacement", "solar altitude"],
        value=default_kinematics_mode,
        description="Kinematics",
    )

    proj_mode_dd = w.Dropdown(
        options=["limb", "custom"],
        value=default_projection_mode,
        description="Projection",
    )
    proj_ang_in = w.FloatText(value=float(default_projection_angle_deg), description="Angle(deg)")

    annotate_dd = w.Dropdown(
        options=["time", "freq", "index"],
        value="time",
        description="Annotate",
    )

    run_btn = w.Button(description="Run Step-6", button_style="success")
    progress = w.IntProgress(
        value=0,
        min=0,
        max=1,
        step=1,
        description="Progress",
        bar_style="",
        orientation="horizontal",
        layout=w.Layout(width="700px"),
    )
    status_html = w.HTML(value="<b>Status:</b> Idle")
    out = w.Output()

    avail_map = {}
    selected_map = {}

    # --------------------------------------------------
    # refresh logic
    # --------------------------------------------------
    def _set_root_from_source(_=None):
        if source_dd.value == "Step5B corrected":
            if step5b_root is not None:
                root_in.value = str(step5b_root)
        else:
            root_in.value = str(step4_root)

    def _update_proj_visibility(*args):
        use_proj = (kin_mode_dd.value == "projected displacement")
        proj_mode_dd.layout.display = None if use_proj else "none"
        proj_ang_in.layout.display = None if use_proj else "none"

    def _update_aia_visibility(*args):
        use_aia = (bg_mode_dd.value == "AIA")
        aia_wave_in.layout.display = None if use_aia else "none"
        aia_time_mode_dd.layout.display = None if use_aia else "none"
        aia_time_manual_in.layout.display = None if (use_aia and aia_time_mode_dd.value == "manual") else "none"
        jsoc_email_in.layout.display = None if use_aia else "none"

    def _refresh(_=None):
        nonlocal avail_map, selected_map
        root_path = Path(root_in.value.strip())
        with out:
            clear_output()

            if not root_path.exists():
                print("Root does not exist:", root_path)
                avail_sel.options = []
                return

            files = _scan_all_image_fits(root_path)
            avail_map = {}
            for f in files:
                rel = str(Path(f).relative_to(root_path))
                avail_map[rel] = f

            avail_opts = sorted(avail_map.keys())
            avail_sel.options = avail_opts

            selected_map = {k: v for k, v in selected_map.items() if Path(v).exists()}
            selected_sel.options = sorted(selected_map.keys())

            print("Source:", source_dd.value)
            print("Root  :", root_path)
            print("Found :", len(files), "candidate FITS")
            if len(files) > 0:
                print("First :", avail_opts[0])
                print("Last  :", avail_opts[-1])

        progress.value = 0
        progress.max = 1
        progress.bar_style = ""
        status_html.value = "<b>Status:</b> Ready"

    def _add_selected(_=None):
        nonlocal selected_map
        for rel in avail_sel.value:
            if rel in avail_map:
                selected_map[rel] = avail_map[rel]
        selected_sel.options = sorted(selected_map.keys())

    def _remove_selected(_=None):
        nonlocal selected_map
        for rel in list(selected_sel.value):
            if rel in selected_map:
                del selected_map[rel]
        selected_sel.options = sorted(selected_map.keys())

    def _add_all(_=None):
        nonlocal selected_map
        for rel, absf in avail_map.items():
            selected_map[rel] = absf
        selected_sel.options = sorted(selected_map.keys())

    def _clear_all(_=None):
        nonlocal selected_map
        selected_map = {}
        selected_sel.options = []

    # --------------------------------------------------
    # run logic
    # --------------------------------------------------
    def _on_run(_):
        with out:
            clear_output()

            root_path = Path(root_in.value.strip())
            run_tag = user_tag_in.value.strip() or "t3beam_run"

            roi = (
                float(xmin_in.value),
                float(xmax_in.value),
                float(ymin_in.value),
                float(ymax_in.value),
            )

            fov = (
                float(fov_xmin_in.value),
                float(fov_xmax_in.value),
                float(fov_ymin_in.value),
                float(fov_ymax_in.value),
            )

            thresh_frac = float(thresh_in.value)
            min_points = int(minpts_in.value)

            bg_mode = bg_mode_dd.value
            aia_wave = int(aia_wave_in.value)
            aia_time_mode = aia_time_mode_dd.value
            aia_time_manual = aia_time_manual_in.value.strip()
            jsoc_email = jsoc_email_in.value.strip()

            kin_mode = kin_mode_dd.value
            proj_mode = proj_mode_dd.value
            proj_ang = float(proj_ang_in.value)
            annotate_by = annotate_dd.value

            files = [selected_map[k] for k in selected_sel.options]
            if len(files) == 0:
                print("No selected FITS.")
                status_html.value = "<b>Status:</b> No selected FITS"
                progress.bar_style = "warning"
                return

            run_out = _ensure_dir(out_root / run_tag)
            csv_path = run_out / "propagation_results.csv"
            jsonl_path = run_out / "propagation_results.jsonl"
            overplot_png = run_out / "overplot_centroids_contours.png"
            kin_png = run_out / "distance_time.png"
            log_path = run_out / "step6.log"
            summary_path = run_out / "summary.txt"

            if overwrite_cb.value:
                for p in [csv_path, jsonl_path, overplot_png, kin_png, summary_path, log_path]:
                    if p.exists():
                        p.unlink()

            progress.max = max(1, len(files))
            progress.value = 0
            progress.bar_style = "info"
            status_html.value = f"<b>Status:</b> Running Step-6 on {len(files)} FITS..."

            print("Source      :", source_dd.value)
            print("Root        :", root_path)
            print("Run tag     :", run_tag)
            print("Out         :", run_out)
            print("N FITS      :", len(files))
            print("ROI         :", roi)
            print("FOV         :", fov)
            print("Background  :", bg_mode)
            if bg_mode == "AIA":
                print("AIA Å       :", aia_wave)
                print("AIA time    :", aia_time_mode, aia_time_manual if aia_time_mode == "manual" else "")
            print("Kinematics  :", kin_mode)
            if kin_mode == "projected displacement":
                print("Projection  :", proj_mode, "angle =", proj_ang)
            print("----")

            results = []
            for i, fpath in enumerate(files, start=1):
                status_html.value = f"<b>Status:</b> Processing {i}/{len(files)} : {os.path.basename(fpath)}"
                try:
                    rec = _measure_one_fits(
                        fits_path=fpath,
                        roi=roi,
                        thresh_frac=thresh_frac,
                        min_points=min_points,
                    )
                    results.append(rec)
                    with open(log_path, "a") as f:
                        f.write(f"OK   {fpath}\n")
                    print("[OK]", rec["sb_tag"], rec["fname"], rec["obstime_isot"])
                except Exception as e:
                    with open(log_path, "a") as f:
                        f.write(f"FAIL {fpath} : {e}\n")
                    print("[FAIL]", os.path.basename(fpath), ":", e)
                progress.value = i

            results = sorted(results, key=lambda r: r["obstime"].unix if "obstime" in r else 1e99)

            aia_map = None
            aia_file = None
            if bg_mode == "AIA":
                try:
                    ref_time = _get_aia_reference_time(
                        results,
                        mode=aia_time_mode,
                        manual_time_str=aia_time_manual,
                    )

                    # save email for later reruns in the same kernel
                    if jsoc_email:
                        os.environ["JSOC_EMAIL"] = jsoc_email

                    # 1) try cache first
                    aia_map, aia_file = _find_cached_aia_map(
                        ref_time=ref_time,
                        wavelength=aia_wave,
                        search_dirs=[run_out, out_root],
                        max_dt=30*u.min,
                    )

                    if aia_map is not None:
                        print("AIA cached  :", aia_file)
                    else:
                        # 2) download only if not cached
                        status_html.value = "<b>Status:</b> Downloading AIA cutout..."
                        aia_map, aia_file = _download_aia_cutout_map(
                            ref_time=ref_time,
                            fov=fov,
                            wavelength=aia_wave,
                            jsoc_email=jsoc_email,
                            out_dir=run_out,
                        )
                        print("AIA file    :", aia_file)

                except Exception as e:
                    print("[WARN] AIA background download failed:", e)
                    print("[WARN] Falling back to blank background.")
                    bg_mode = "blank"

            results, fit_meta = _add_kinematics(
                results,
                kinematics_mode=kin_mode,
                projection_mode=proj_mode,
                angle_deg=proj_ang,
            )

            with open(csv_path, "w") as f:
                f.write(
                    "sb_tag,fname,obstime_isot,freq_mhz,ok,"
                    "cen_Tx_as,cen_Ty_as,fwhm_major_as,fwhm_minor_as,"
                    "rsun_obs_as,r_as,altitude_as,altitude_km,altitude_Mm,"
                    "dt_s,proj_coord_as,proj_coord_km,proj_coord_Mm,"
                    "pairwise_speed_km_s,pairwise_speed_c,"
                    "pairwise_altitude_speed_km_s,pairwise_altitude_speed_c,"
                    "fit_speed_km_s,fit_speed_c,altitude_fit_speed_km_s,altitude_fit_speed_c\n"
                )
                for r in results:
                    f.write(
                        f"{r.get('sb_tag','')},{r.get('fname','')},{r.get('obstime_isot','')},{r.get('freq_mhz',np.nan)},{int(r.get('ok',False))},"
                        f"{r.get('cen_Tx_as',np.nan)},{r.get('cen_Ty_as',np.nan)},{r.get('fwhm_major_as',np.nan)},{r.get('fwhm_minor_as',np.nan)},"
                        f"{r.get('rsun_obs_as',np.nan)},{r.get('r_as',np.nan)},{r.get('altitude_as',np.nan)},{r.get('altitude_km',np.nan)},{r.get('altitude_Mm',np.nan)},"
                        f"{r.get('dt_s',np.nan)},{r.get('proj_coord_as',np.nan)},{r.get('proj_coord_km',np.nan)},{r.get('proj_coord_Mm',np.nan)},"
                        f"{r.get('pairwise_speed_km_s',np.nan)},{r.get('pairwise_speed_c',np.nan)},"
                        f"{r.get('pairwise_altitude_speed_km_s',np.nan)},{r.get('pairwise_altitude_speed_c',np.nan)},"
                        f"{r.get('fit_speed_km_s',np.nan)},{r.get('fit_speed_c',np.nan)},{r.get('altitude_fit_speed_km_s',np.nan)},{r.get('altitude_fit_speed_c',np.nan)}\n"
                    )

            with open(jsonl_path, "w") as f:
                for r in results:
                    rr = dict(r)
                    if "obstime" in rr:
                        rr["obstime"] = rr["obstime"].isot
                    rr = _jsonify_obj(rr)
                    f.write(json.dumps(rr) + "\n")

            _plot_overplot(
                results,
                fov,
                overplot_png,
                annotate_by=annotate_by,
                background_mode=bg_mode,
                aia_map=aia_map,
            )
            _plot_kinematics(results, fit_meta, kin_png)

            good = [r for r in results if r.get("ok", False)]
            with open(summary_path, "w") as f:
                f.write(f"run_tag: {run_tag}\n")
                f.write(f"source: {source_dd.value}\n")
                f.write(f"root: {root_path}\n")
                f.write(f"n_input: {len(files)}\n")
                f.write(f"n_ok: {len(good)}\n")
                f.write(f"roi: {roi}\n")
                f.write(f"fov: {fov}\n")
                f.write(f"background_mode: {bg_mode}\n")
                if aia_file is not None:
                    f.write(f"aia_file: {aia_file}\n")
                f.write(f"kinematics_mode: {kin_mode}\n")
                if kin_mode == "projected displacement":
                    f.write(f"projection_mode: {proj_mode}\n")
                    f.write(f"projection_angle_deg: {proj_ang}\n")
                if fit_meta is not None:
                    f.write(f"fit_speed_km_s: {fit_meta.get('fit_speed_km_s', np.nan)}\n")
                    f.write(f"fit_speed_c: {fit_meta.get('fit_speed_c', np.nan)}\n")
                    f.write(f"fit_speed_err_km_s: {fit_meta.get('fit_speed_err_km_s', np.nan)}\n")
                    f.write(f"fit_speed_err_c: {fit_meta.get('fit_speed_err_c', np.nan)}\n")
                    f.write(f"pairwise_mean_km_s: {fit_meta.get('pairwise_mean_km_s', np.nan)}\n")
                    f.write(f"pairwise_std_km_s: {fit_meta.get('pairwise_std_km_s', np.nan)}\n")
                    f.write(f"pairwise_mean_c: {fit_meta.get('pairwise_mean_c', np.nan)}\n")
                    f.write(f"pairwise_std_c: {fit_meta.get('pairwise_std_c', np.nan)}\n")
                    f.write(f"n_unique_times: {fit_meta.get('n_unique_times', np.nan)}\n")

            progress.value = progress.max
            progress.bar_style = "success"
            status_html.value = "<b>Status:</b> Done"

            print("----")
            print("Done.")
            print("CSV     :", csv_path)
            print("JSONL   :", jsonl_path)
            print("Overplot:", overplot_png)
            print("Kine fig:", kin_png)
            print("Summary :", summary_path)
            print("Log     :", log_path)

            if fit_meta is not None:
                if np.isfinite(fit_meta.get("fit_speed_km_s", np.nan)):
                    msg = f'Linear-fit speed: {fit_meta["fit_speed_km_s"]:.1f} km/s ({fit_meta["fit_speed_c"]:.4f} c)'
                    if np.isfinite(fit_meta.get("fit_speed_err_km_s", np.nan)):
                        msg += f' ± {fit_meta["fit_speed_err_km_s"]:.1f} km/s'
                    print(msg)

                if np.isfinite(fit_meta.get("pairwise_mean_km_s", np.nan)):
                    msg = (
                        f'Pairwise speed mean ± std: '
                        f'{fit_meta["pairwise_mean_km_s"]:.1f} ± {fit_meta["pairwise_std_km_s"]:.1f} km/s '
                        f'({fit_meta["pairwise_mean_c"]:.4f} ± {fit_meta["pairwise_std_c"]:.4f} c)'
                    )
                    print(msg)

            if overplot_png.exists():
                print("\n[Preview] Overplot figure")
                display(IPyImage(filename=str(overplot_png)))
            if kin_png.exists():
                print("\n[Preview] Kinematics figure")
                display(IPyImage(filename=str(kin_png)))

    # --------------------------------------------------
    # wiring
    # --------------------------------------------------
    refresh_btn.on_click(_refresh)
    add_btn.on_click(_add_selected)
    remove_btn.on_click(_remove_selected)
    add_all_btn.on_click(_add_all)
    clear_btn.on_click(_clear_all)
    run_btn.on_click(_on_run)
    source_dd.observe(_set_root_from_source, names="value")
    kin_mode_dd.observe(_update_proj_visibility, names="value")
    bg_mode_dd.observe(_update_aia_visibility, names="value")
    aia_time_mode_dd.observe(_update_aia_visibility, names="value")

    _refresh()
    _update_proj_visibility()
    _update_aia_visibility()

    display(
        w.VBox([
            w.HBox([source_dd, refresh_btn, overwrite_cb]),
            root_in,
            w.HBox([
                avail_sel,
                w.VBox([add_btn, remove_btn, add_all_btn, clear_btn]),
                selected_sel,
            ]),
            w.HBox([user_tag_in]),
            w.HBox([xmin_in, xmax_in, ymin_in, ymax_in]),
            w.HBox([fov_xmin_in, fov_xmax_in, fov_ymin_in, fov_ymax_in]),
            w.HBox([thresh_in, minpts_in]),
            w.HBox([bg_mode_dd, aia_wave_in, aia_time_mode_dd]),
            w.HBox([aia_time_manual_in, jsoc_email_in]),
            w.HBox([kin_mode_dd, proj_mode_dd, proj_ang_in]),
            w.HBox([annotate_dd]),
            run_btn,
            progress,
            status_html,
            out,
        ])
    )