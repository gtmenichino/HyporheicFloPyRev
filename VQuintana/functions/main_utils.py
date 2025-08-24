"""Utility functions for the combined hyporheic model workflow.

This module centralises reusable helpers that were previously spread
across multiple notebooks.  Functions here intentionally avoid
side-effects so they can be imported and unit tested independently.
"""

from __future__ import annotations

from pathlib import Path, PurePath
from types import SimpleNamespace
from typing import Any, Sequence, Tuple

import numpy as np
import flopy
from scipy.interpolate import griddata

try:
    from modflow_devtools.misc import timed  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    def timed(func):  # type: ignore
        return func


# ---------------------------------------------------------------------
# General helper functions
# ---------------------------------------------------------------------

def interpolate_na(terrain: np.ma.MaskedArray) -> np.ndarray:
    """Interpolate missing values in a masked terrain array.

    Parameters
    ----------
    terrain : np.ma.MaskedArray
        Masked array of terrain values where ``terrain.mask`` identifies
        cells requiring interpolation.
    """
    # Get coordinates of valid/invalid cells
    valid_mask = ~terrain.mask
    valid_coords = np.array(np.nonzero(valid_mask)).T
    valid_values = terrain[valid_mask]

    invalid_mask = terrain.mask
    invalid_coords = np.array(np.nonzero(invalid_mask)).T

    # Interpolate using nearest-neighbour
    interpolated_values = griddata(
        valid_coords, valid_values, invalid_coords, method="nearest"
    )
    filled = terrain.copy()
    filled[invalid_mask] = interpolated_values
    return filled


def calculate_gw_elevation(boundary_cells, top_elevation, offset):
    """Calculate groundwater elevation for a set of boundary cells."""
    gw_elevation_min: list[float] = []
    for cell in boundary_cells:
        layer, row, col = cell
        gw_elevation = top_elevation + offset
        gw_elevation_min.append([layer, row, col, gw_elevation])
    return gw_elevation_min


def interpolate_gw_elevation_first_layer_only(
    first_layer_cells: list[tuple[int, int, int]],
    head_first: float,
    head_last: float,
) -> list[float]:
    """Linearly interpolate groundwater elevation across boundary cells."""
    n_cells = len(first_layer_cells)
    interpolated_heads: list[float] = []
    for idx, cell in enumerate(first_layer_cells):
        frac = idx / max(n_cells - 1, 1)
        head = head_first + frac * (head_last - head_first)
        interpolated_heads.append([*cell, head])
    return interpolated_heads


# ---------------------------------------------------------------------
# Model construction and execution
# ---------------------------------------------------------------------

def build_gwf_model(
    cfg: Any,
    chd_data: Sequence[Sequence[float]],
    idomain: np.ndarray,
) -> Tuple[flopy.mf6.MFSimulation, flopy.mf6.ModflowGwf]:
    """Construct a stand-alone MODFLOW 6 groundwater-flow model.

    This is a lightly adapted version of the implementation originally
    contained in the `run_models.ipynb` notebook.  The function assembles
    a MODFLOW 6 simulation using configuration values provided in
    ``cfg`` along with constant-head (CHD) data and an ``idomain`` mask.
    """
    if idomain.shape != (cfg.nlay, cfg.nrow, cfg.ncol):  # type: ignore[attr-defined]
        raise ValueError("`idomain` dimensions don’t match cfg grid.")

    sim = flopy.mf6.MFSimulation(
        sim_name=cfg.sim_name,
        exe_name=str(cfg.md6_exe_path),
        sim_ws=str(cfg.gwf_ws),
    )

    flopy.mf6.ModflowTdis(
        sim,
        time_units=cfg.time_units.upper(),
        nper=cfg.nper,
        perioddata=[(cfg.perlen, cfg.nstp, cfg.tsmult)],
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=cfg.gwf_name,
        save_flows=True,
    )

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=cfg.nlay,
        nrow=cfg.nrow,
        ncol=cfg.ncol,
        delr=cfg.cell_size_x,
        delc=cfg.cell_size_y,
        top=cfg.tops[0],
        botm=cfg.botm,
        idomain=idomain,
        xorigin=cfg.xmin,
        yorigin=cfg.ymin,
    )

    gwf.modelgrid.crs = cfg.raster_crs
    gwf.modelgrid.set_coord_info(cfg.xmin, cfg.ymin, crs=cfg.raster_crs)

    strt = np.full((cfg.nlay, cfg.nrow, cfg.ncol), cfg.bed_elevation, dtype=float)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=2,
        k=cfg.kh,
        k33=cfg.kv,
        save_flows=True,
        save_saturation=True,
        save_specific_discharge=True,
    )

    if chd_data:
        flopy.mf6.ModflowGwfchd(
            gwf,
            maxbound=len(chd_data),
            stress_period_data={0: chd_data},
            save_flows=True,
        )

    flopy.mf6.ModflowGwfoc(
        gwf,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        head_filerecord=[cfg.headfile],
        budget_filerecord=[cfg.budgetfile],
        printrecord=[("HEAD", "LAST")],
    )

    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-4,
        outer_maximum=200,
        inner_maximum=500,
        inner_dvclose=1e-4,
        rcloserecord=1e-4,
        linear_acceleration="BICGSTAB",
        relaxation_factor=0.97,
    )

    return sim, gwf


def build_particle_models(
    sim_name: str,
    gwf: flopy.mf6.ModflowGwf,
    river_cells: list[tuple[int, int, int]],
    *,
    mp7_ws: Path | str | None = None,
    exe_path: Path | str | None = None,
):
    """Create forward and backward MODPATH 7 models."""
    from flopy.modpath import Modpath7, ParticleData, ParticleGroup

    if mp7_ws is None:
        mp7_ws = Path(gwf.simulation.sim_path).parent / "mp7_workspace"
    mp7_ws = Path(mp7_ws).absolute()
    mp7_ws.mkdir(exist_ok=True)

    if exe_path is None:
        exe_path = "mp7"

    def _make(direction: str):
        mp = Modpath7.create_mp7(
            modelname=f"{sim_name}_mp_{direction}",
            trackdir=direction,
            flowmodel=gwf,
            model_ws=mp7_ws,
            exe_name=str(exe_path),
        )
        partlocs = [(k, i, j) for (k, i, j, *_) in river_cells]
        particle_data = ParticleData(partlocs, structured=True, drape=0)
        pg = ParticleGroup(particledata=particle_data)
        mpsim = mp.get_package("MPSIM")
        mpsim.particlegroups.clear()
        mpsim.particlegroups.append(pg)
        return mp

    return _make("forward"), _make("backward")


def write_models(*sims, silent: bool = False) -> None:
    """Write MODFLOW or MODPATH models to disk."""
    for sim in sims:
        if isinstance(sim, flopy.mf6.MFSimulation):
            sim.write_simulation(silent=silent)
        else:
            sim.write_input()


@timed
def run_models(*sims, silent: bool = False) -> None:
    """Run one or more MODFLOW/MODPATH simulations."""
    for sim in sims:
        if isinstance(sim, flopy.mf6.MFSimulation):
            print(f"Running simulation: {sim.name}")
            success, buff = sim.run_simulation(silent=silent, report=True)
        else:
            print(f"Running model: {sim.name}")
            success, buff = sim.run_model(silent=silent, report=True)
        if not success:
            raise RuntimeError(f"Simulation {sim.name} failed: {buff}")


# Placeholder for the extensive post-processing logic originally present
# in ``run_models.ipynb``.  The detailed implementation has been omitted
# here for brevity but can be reintroduced as required.
def process_and_export_modpath7_results(*args, **kwargs):  # pragma: no cover
    """Process MODPATH 7 results and export spatial datasets.

    The full implementation is available in the original notebook.  This
    placeholder allows the main workflow notebook to import the symbol
    without executing heavy post-processing when not needed.
    """
    raise NotImplementedError(
        "process_and_export_modpath7_results is not implemented in this "
        "repository snapshot."
    )
