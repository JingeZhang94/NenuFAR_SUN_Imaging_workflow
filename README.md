[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18848880.svg)](https://doi.org/10.5281/zenodo.18848880)

# NenuFAR_SUN_Imaging_workflow

Reproducible NenuFAR solar interferometric imaging workflow **(Step1–4)** plus **Step5 ionospheric-offset correction** and **centroid / quicklook** tools.

This repository provides a lightweight, notebook-driven UI (via `ipywidgets`) to:
- select NenuFAR solar sub-bands (SBs) and associated CASA calibrator MS products,
- run standard imaging steps (**DP3 + WSClean**),
- generate quicklooks / movies,
- derive a quiet-Sun ionospheric offset solution (**Step5A**),
- apply the solution to selected/all FITS (**Step5B**, producing corrected FITS + before/after quicklooks),
- measure source centroids with a 2D Gaussian fit in an ROI (**Step5C**, saving tables/records + optional movie).

> **Scope:** This is a *workflow + UI* for producing solar imaging products and diagnostic quicklooks. The emphasis is on reproducibility and practical usage on the NenuFAR computing environment.

---

## Important requirements (read first)

### 1) Where to run / access policy (Nançay servers)
This workflow is intended to run on the **Nançay computing environment** (e.g., `nancep`), where NenuFAR data paths and the standard processing toolchain are available.

Access to the relevant data and compute environment may require approval from the **NenuFAR KP11 (Solar) team**.  
For permission requests, please contact the current PI: **Carine Briand** (carine.briand@obspm.fr).

### 2) Two external files are required (not included in this repo)

**(a) Container image**  
`linc_latest.sif`  
Contains the runtime environment used for DP3/WSClean execution.

**(b) Calibration database**  
`CasA.sourcedb`  
Used by the calibration steps.

If you do not have access to these files (or need verified paths), **please contact the author**.

---

## Requirements

### Python environment
Tested with Python 3.8+ (recommend 3.12.3 typical for SunPy/Astropy stacks). Required Python packages:
- `numpy`, `matplotlib`
- `astropy`
- `sunpy`
- `ipywidgets`
- `scipy` (for 2D Gaussian fitting)

Optional:
- `ffmpeg` (recommended) for MP4 movie generation; otherwise GIF fallback may be used.

### External tools (HPC / cluster)
For Step1–4 you need the standard imaging toolchain available on your system, e.g.
- `DP3`
- `wsclean`
- container runtime (e.g. `apptainer` / `singularity`)

> Step5 (A/B/C) operates on **FITS images** produced by Step4 and does not require DP3/WSClean.

---

## Repository layout

Main entry points:
- `nenufar_ui.py` — notebook UI functions (Step1–5)
- `workflow_v020.ipynb` — example notebook showing how to run the workflow
- `nenufar_sb_scan.py` — helper for scanning candidate SBs / associations

---

## Quick start

### 1) Clone the repository
```bash
git clone https://github.com/JingeZhang94/NenuFAR_SUN_Imaging_workflow.git
cd NenuFAR_SUN_Imaging_workflow
```

---

### 2) Open the example notebook

On the Nançay / `nancep` environment, open the example notebook:

- `workflow_v020.ipynb`

This notebook is the recommended entry point. You will run the workflow by executing a small number of cells (Step1 → Step5).  
Each step launches an interactive UI (`ipywidgets`): you set the paths/parameters in the widgets, then click **Run**.

**Typical usage pattern**
1. Open `workflow_v020.ipynb` (JupyterLab or VSCode Jupyter both work).
2. Run cells from top to bottom.
3. In each step UI:
   - choose the SB(s) / FITS,
   - adjust parameters if needed,
   - click **Run**,
   - check the printed log + generated outputs under the corresponding `step*_outputs_YYYYMMDD/` folder.

If you only want to work on ionospheric correction / centroid tools, you can start from **Step5** directly as long as you already have Step4 FITS images available.

### 3) Run Step1–Step4 (imaging + quicklook)

All Step1–Step4 operations are driven from the example notebook `workflow_v020.ipynb`.
You typically run the notebook **top → bottom**, and each step opens an interactive UI (`ipywidgets`).

#### 3.1 (Optional) Scan & select SB candidates (Nançay / nancep only)
If you want to scan available NenuFAR SBs and quickly locate relevant sub-bands:

- Open the Step0/scan cell in `workflow_v020.ipynb`
- Set the base path(s) (e.g. NenuFAR data root) and work directory
- Click **Run** in the UI
- The table output can be used to identify SBs to process

> Note: this scan utility assumes the Nançay / `nancep` computing environment and access permissions to NenuFAR data paths.

#### 3.2 Step1 — Prepare MS (CASA pre/post) per SB
- Run the Step1 cell (calls `nenufar_ui.run_step1_ui(...)`)
- In the UI:
  - choose SB(s)
  - choose calibration mode (dropdown / pre / post)
  - set output root (recommended: `.../step1_outputs_YYYYMMDD`)
  - click **Run Step-1**
- Outputs:
  - per-SB folder under `step1_outputs_YYYYMMDD/`
  - log printed in the notebook + written to the output folder

#### 3.3 Step2 — Time ROI cut (make ROI MS)
- Run the Step2 cell (calls `nenufar_ui.run_step2_zoom_ui(...)`)
- In the UI:
  - set Start / End time range
  - select SB(s)
  - click **Run Step-2 (ROI)**
- Outputs:
  - ROI MS under `step2_outputs_YYYYMMDD/ROI/SBxxx/`
  - log printed + saved

#### 3.4 Step3 — Calibration (DP3)
- Run the Step3 cell (calls `nenufar_ui.run_step3_calib_ui(...)`)
- In the UI:
  - select SB(s)
  - click **Run Step-3 (Calib)**
- Outputs:
  - per-SB calibration products under `step3_outputs_YYYYMMDD/SBxxx/`
  - parset files + log

#### 3.5 Step4 — Imaging (WSClean) + quicklooks
- Run the Step4 imaging cell (calls `nenufar_ui.run_step4_wsclean_ui(...)`)
- In the UI:
  - set WSClean parameters (image size, scale, weighting, robust, etc.)
  - select SB(s)
  - click **Run Step-4 (WSClean)**
- Outputs:
  - FITS images under `step4_outputs_YYYYMMDD/SBxxx/*-image.fits`

Then generate Step4-style PNG quicklooks / movies:
- Run the Step4 quicklook cell (calls `nenufar_ui.run_step4_quicklook_ui(...)`)
- In the UI:
  - select one or multiple `*-image.fits`
  - set crop size + color limits + contours
  - click **Generate Quicklooks**
- Outputs:
  - `step4_outputs_YYYYMMDD/quicklook/SBxxx/*.png`
  - optional movie (`*.mp4` or `*.gif`) + `quicklook.log`

### 4) Run Step5A / Step5B / Step5C (ionospheric offset correction + centroid tool)

Step5 operates on FITS images (typically produced by Step4).

#### Step5A — solve an ionospheric offset from a Quiet Sun frame
- Run the **Step5A** cell in the notebook.
- In the widget:
  - choose the SB,
  - choose a *Quiet Sun* FITS frame (before the burst onset),
  - set an ROI (arcsec) that isolates the Quiet Sun emission,
  - click **Solve Step5A**.
- Output (a small JSON solution + diagnostic figure) is written under:
  - `step_iocorrect_outputs_YYYYMMDD/SBxxx/step5a_solution.json`

The solution contains the WCS shift to apply (in pixel space on the Step4 map grid).

#### Step5B — apply the Step5A solution to FITS (write corrected FITS + before/after quicklooks)
- Run the **Step5B** cell.
- In the widget:
  - choose the SB,
  - either:
    - select specific FITS files with the mouse, **or**
    - enable **Select ALL image.fits** to batch-process everything,
  - click **Run Step5B (apply)**.
- Outputs are written under:
  - corrected FITS: `step_iocorrect_outputs_YYYYMMDD/SBxxx/corr_fits/*.fits`
  - before/after quicklooks: `step_iocorrect_outputs_YYYYMMDD/SBxxx/quicklook_step5b/`
  - log: `step_iocorrect_outputs_YYYYMMDD/SBxxx/step5b_apply.log`

#### Step5C — centroid extraction (works for raw Step4 FITS or corrected Step5B FITS)
- Run the **Step5C** cell.
- In the widget:
  - choose **Source**: `Step4 raw` or `Step5B corrected`,
  - select FITS file(s) (or select all),
  - define an ROI (arcsec) around the radio source,
  - set `thresh_frac` and `min_points` for the 2D Gaussian fit,
  - click **Run Step5C (centroid)**.
- Outputs are written under:
  - quicklook PNGs: `step_iocentroid_outputs_YYYYMMDD/.../quicklook_centroid/`
  - table (CSV): `centroid_results.csv`
  - records (JSONL): `centroid_results.jsonl`
  - log: `step5c_centroid.log`
  - optional movie (if enabled): `centroid_movie.mp4` (or GIF fallback)

The centroid table/records include (when available): date tag, SB, FITS filename (t-index), observation time, ROI bounds, and centroid coordinates (arcsec / world coordinates).

---

## Acknowledgements

This project has received funding from the European Union's Horizon Europe research and innovation programme under grant agreement No 101134999. This repository reflects only the author's view and the European Commission is not responsible for any use that may be made of the information it contains.


**Cite (all versions):** https://doi.org/10.5281/zenodo.18848880  