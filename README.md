[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18848880.svg)](https://doi.org/10.5281/zenodo.18848880)

# NenuFAR_SUN_Imaging_workflow

**`SUN_Imaging.ipynb`** is the main notebook interface of this workflow.

It provides a **simple and efficient way** to process **NenuFAR solar interferometric imaging data**, including the main imaging workflow and the most commonly used post-processing tools.

It is intended to make the processing chain more accessible and reproducible, including for **non-specialist users**.

---

## Overview

This repository provides a lightweight UI-based workflow (via `ipywidgets`) to:

- select NenuFAR solar sub-bands (SBs) and associated calibrator products,
- prepare and process raw interferometric data,
- run the main imaging chain based on **DP3 + WSClean**,
- generate FITS quicklooks / movies,
- derive and apply WCS offset corrections from Quiet Sun frames,
- measure source centroids and uncertainties,
- perform beam-propagation / source-tracking analysis from selected FITS images.

> **Important:**  
> The **main imaging workflow (Step1–Step4)** is intended to run in the **Nançay computing environment** (in particular on **`nancep`**), because the raw NenuFAR interferometric **MS data** are stored there and are typically too large to be practically downloaded and reprocessed elsewhere in a conventional way.  
> The **post-processing tools** can also be used separately on selected FITS products.

---

## Main workflow vs post-processing tools

### Main imaging workflow

The core imaging workflow consists of four steps:

- **Step1 — Prepare MS**  
  Prepare calibrator and solar MS products for further processing.

- **Step2 — ROI extraction**  
  Extract the solar time window / ROI from the prepared solar MS.

- **Step3 — Calibration transfer**  
  Derive gain solutions from the calibrator data and apply them to the solar ROI MS.

- **Step4 — Imaging**  
  Run WSClean to generate FITS image products, and optionally generate quicklook PNGs / movies.

These tools are provided in:

- `nenufar_ui.py`

### Post-processing tools

The workflow also includes a set of post-processing tools that operate on FITS products:

- **Post-processing tool A — WCS offset solve**  
  Derive a WCS offset solution from a Quiet Sun frame.

- **Post-processing tool B — WCS correction apply**  
  Apply the saved WCS offset solution to selected FITS frames.

- **Post-processing tool C — Centroid measurement**  
  Measure source centroids and uncertainties from selected FITS frames.

- **Post-processing tool D — Beam propagation / projected source-height analysis**  
  Use measured source positions to derive propagation diagnostics such as projected displacement or solar altitude.

These tools are provided in:

- `post_analysis_tools.py`

---

## Repository layout

Main files:

- `SUN_Imaging.ipynb` — example notebook and recommended entry point
- `nenufar_ui.py` — main imaging workflow tools (**Step1–Step4**)
- `post_analysis_tools.py` — post-processing tools
- `nenufar_sb_scan.py` — helper for scanning candidate sub-bands and associated products
- `README.md`
- `CITATION.cff`
- `LICENSE`

Calibrator sky models:

- `CasA.skymodel` — sky model for **Cassiopeia A**
- `CygA.skymodel` — sky model for **Cygnus A**
- `VirA.skymodel` — sky model for **Virgo A**

Local processing environment:

- `linc_latest.sif` — local container image used for the main imaging workflow, including **DP3** and **WSClean**.  
  This file is **not included** in the GitHub repository and should be obtained separately. See **Important requirements (read first)** below.
---

## Important requirements (read first)

The full imaging workflow (**Step1–Step4**) is intended to run in the **Nançay computing environment**, in particular on **`nancep`**.

The main reason is practical: the raw NenuFAR interferometric **MS data** are stored there and are typically very large, so they are **not meant to be downloaded and reprocessed elsewhere in a conventional way**. The workflow is therefore designed to run **close to the data**.

In the current workflow setup:

- the raw NenuFAR **MS data** are available at Nançay
- the required calibrator sky models are provided in this repository
- the main imaging software environment is provided through the local container **`linc_latest.sif`**, which includes **DP3** and **WSClean**

Access to the relevant data and computing environment may still require approval from the **NenuFAR KP11 (Solar) team**.

For access / permission questions, please contact the current PI:

- **Carine Briand** — `carine.briand@obspm.fr`

For users outside Nançay, the post-processing tools can still be used separately on selected FITS products.

### Container for the main imaging workflow

The full **Step1–Step4** imaging chain uses the container:

- `linc_latest.sif`

This container includes the main processing tools, in particular **DP3** and **WSClean**.

Because of its size, this file is **not distributed through the GitHub repository**.

A local copy is currently available on **`nancep`** at:

- `/data/jzhang/linc_latest.sif`

Users on `nancep` can copy this file into their own repository / working directory if needed.

The container source referenced here follows the LOFAR data-processing tutorial materials by **Harish Vedantham**, and the related tutorial resource is gratefully acknowledged here：

- https://stellar-h2020.eu/index.php/2023/05/23/introduction-to-lofar-data-processing-tutorial/

Example command:

```bash
singularity pull docker://astronrd/linc
```

---

## Installation requirements

To run **`SUN_Imaging.ipynb`**, a standard Python notebook environment is needed.

Tested with:

- **Python 3.12.3**
- **JupyterLab** / notebook
- `ipywidgets` enabled

### Required Python packages

Install the following packages:

```bash
pip install jupyterlab notebook ipywidgets ipympl
pip install numpy pandas scipy matplotlib imageio
pip install astropy sunpy reproject
pip install python-casacore
```

### Example local environment

A simple example using `venv`:

```bash
python -m venv nenufar_sun_env
source nenufar_sun_env/bin/activate
pip install --upgrade pip
pip install jupyterlab notebook ipywidgets ipympl
pip install numpy pandas scipy matplotlib imageio
pip install astropy sunpy reproject
pip install python-casacore
```

Then launch:

```bash
jupyter lab
```

and open:

- `SUN_Imaging.ipynb`

> **Note:**  
> The notebook interface itself can be set up in a standard Python environment.  
> However, the full **Step1–Step4** imaging workflow is still intended to run in the **Nançay / nancep** environment, where the raw NenuFAR **MS data** are stored.

---

## Quick start

### 1) Log in to the Nançay environment

Make sure you are logged in to the **`nancep`** server and working in your own workspace.

This is important because the raw NenuFAR interferometric **MS data** are stored there, and the main imaging workflow is designed to run close to the data.

### 2) Clone the repository

```bash
git clone https://github.com/JingeZhang94/NenuFAR_SUN_Imaging_workflow.git
cd NenuFAR_SUN_Imaging_workflow
```

### 3) Prepare the local workflow environment

### 3) Prepare the local workflow environment

1. Create and activate a Python environment, then install the required packages listed above.
2. Make sure the required container `linc_latest.sif` is available locally, as noted above (see **Important requirements (read first)**).

### 4) Open the notebook

Launch JupyterLab:

```bash
jupyter lab
```

and open:

- `SUN_Imaging.ipynb`

### 5) Set your output paths before running

Before running the workflow, make sure the second configuration cell in the notebook is edited correctly for your own setup, especially the output paths.

This is important because the notebook writes workflow products to the paths defined there.


### 6)Run the notebook step by step

Run the notebook cells one by one, from top to bottom.

Typical usage:

1. configure the paths in the notebook
2. load and inspect the selected data
3. launch each UI from the notebook
4. inspect outputs in the corresponding output folders

---

## Credit, ownership, and responsible use

This workflow was independently developed by **Dr. Jinge Zhang** during his postdoctoral appointment at the **Observatoire de Paris** in 2024-2026.

The code and associated intellectual property belong to the **Observatoire de Paris**, within the framework of the relevant institutional and project context. In addition, a version of this workflow is intended to be distributed and maintained within the **SOLER / Observatoire GitLab infrastructure** for institutional use, long-term maintenance, and future updates.

If this workflow contributes to scientific results, derived data products, or future software developments, users are expected to:

- **cite this workflow and its software release**
- acknowledge its contribution appropriately in related publications, presentations, and derived projects

This repository is shared to support scientific use, reproducibility, and collaboration. These principles are also fully aligned with the spirit of the **SOLER project** and with the original purpose for which this workflow was developed. Its use should respect normal standards of academic attribution, software credit, and intellectual-property respect.

Substantial reuse, repackaging, reshaping, or extension of this workflow without proper attribution is strongly discouraged, especially when this involves AI-assisted restructuring or redevelopment based on the existing concepts, design, or implementation of the tool.

---

## Citation

If you use this workflow in scientific work, please cite the repository / software release using the information provided in:

- `CITATION.cff`

Please also cite the Zenodo record:

- **Cite (all versions):** https://doi.org/10.5281/zenodo.18848880

---

## License

See:

- `LICENSE`

---

## Acknowledgements

This project has received funding from the European Union's Horizon Europe research and innovation programme under grant agreement No 101134999. This repository reflects only the author's view and the European Commission is not responsible for any use that may be made of the information it contains.

Support from **Carine Briand** (NenuFAR KP11 PI), as well as discussions and suggestions from **Alan Loh**, **Julien Girard**, and **Shilpi Bhunia**, are gratefully acknowledged.
