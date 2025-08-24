# Hyporheic FloPy — Rapid Hyporheic‐Flow Modeling Toolkit

Welcome to **Hyporheic FloPy**, a toolkit for **rapidly modeling hyporheic‑zone hydraulics in fluvial systems**.  The project bundles a comment‑preserving YAML configuration, a click‑and‑explore GUI, and ready‑to‑run Jupyter notebooks so you can go from raw spatial data to calibrated MODFLOW 6 results in minutes rather than days.

> **Why model the hyporheic zone?**  Exchange between surface water and shallow groundwater drives temperature buffering, nutrient cycling, and ecological refuge in streams.  Quantifying that exchange helps engineers and scientists design more resilient, nature‑based restoration measures.

---

## Key capabilities

| Feature                                   | What it does                                                                                               | Where to start                              |
| ----------------------------------------- | ---------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| **Automated workspace builder**           | Converts rasters & shapefiles into a fully meshed MODFLOW 6 model (including MODPATH 7 particle tracking). | `__main__.py` or the *Quick‑start notebook* |
| **Comment‑preserving YAML inputs**        | Single‑file config with sensible defaults that you can edit in‑app or in your editor of choice.            | `inputs.yaml`                               |
| **PyQt GUI**                              | Inspect layers, toggle boundaries, zoom, and launch model runs without leaving VS Code.                    | `app_gui.py`                                |
| **Notebooks for exploration & reporting** | Re‑usable notebooks show every processing step and generate a multi‑page PDF summary.                      | `notebooks/`                                |
| **Modular helpers**                       | `functions/*` packages split raster, vector, and MODFLOW logic for straightforward unit testing.           | Dive into `functions/`                      |

---

## Quick start

### 1 . Clone & install

```bash
# clone your fork
$ git clone https://github.com/<you>/HyporheicFloPy.git
$ cd HyporheicFloPy

# (optional) set up an isolated env
$ python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate

# install core deps – everything else is pulled on‑demand
$ pip install -r requirements.txt
```

> **Heads‑up:**  The first run downloads official USGS MODFLOW 6 and MODPATH 7 binaries into `modflowExe/` (≈ 40 MB).  Subsequent runs reuse the local copies.

### 2 . Launch the GUI *(best for first‑time users)*

```bash
$ python app_gui.py
```

* Open the repo in **Visual Studio Code**.
* Press F5 (or use the *Run ▶* button) to debug `app_gui.py` step‑by‑step.
* Tweak any field in **Setup → inputs.yaml**, hit **Geometry → Load/refresh**, and watch the grid rebuild in real time.

### 3 . Run the notebooks *(transparent, reproducible workflow)*

The `notebooks/` folder contains a linear, bite‑sized workflow:

1. **00\_setup.ipynb** – environment check & helper imports
2. **10\_build\_grid.ipynb** – create the finite‑difference grid
3. **20\_boundary\_conditions.ipynb** – derive river & CHD inputs
4. **30\_run\_models.ipynb** – write, run, and post‑process MF6/MP7

Feel free to copy a notebook and swap in your own rasters or shapefiles.

---

## Repository layout

```
.
├─ app_gui.py            # PyQt GUI entry‑point
├─ __main__.py           # CLI pipeline (python -m HyporheicFloPy)
├─ inputs.yaml           # One‑stop project configuration
├─ notebooks/            # Jupyter notebooks (tutorial & benchmarking)
├─ functions/            # Re‑usable helpers (raster_utils, model_utils …)
├─ modflowExe/           # MODFLOW 6 & MODPATH 7 binaries (auto‑downloaded)
└─ VQuintana/            # Example data bundle (rasters, shapefiles)
```

---

## Citing & contributing

* **Citation:**  If Hyporheic FloPy aided your research, please cite the companion preprint: *Menichino et al., 2025, “Accelerating Hyporheic‑Zone Modeling with FLOPY and Python Workflows”.*

* **License:**  MIT; see `LICENSE` for details.

---