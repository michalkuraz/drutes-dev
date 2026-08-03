"""Schema and factory for ``drutes.conf/global.conf``."""

from __future__ import annotations

from pathlib import Path

from .configfile import ConfigFile
from .parameter import ParameterDefinition as Definition
from .parameter import ParameterType as Type

GLOBAL_DEFINITIONS = (
    Definition("model_type", "Model type", Type.CHOICE, choices=("RE", "REstd", "boussi", "ADE", "ADEnc", "Re_dual", "heat")),
    Definition("dimension", "Problem dimension", Type.CHOICE, choices=("1", "2", "2r", "3")),
    Definition("mesh_generator", "Mesh generator", Type.INTEGER),
    Definition("max_picard_iterations", "Maximum Picard iterations", Type.INTEGER),
    Definition("h_tolerance", "Picard iteration tolerance", Type.FLOAT),
    Definition("time_units", "Time units", help_text="Up to five characters."),
    Definition("dt", "Initial time step", Type.FLOAT),
    Definition("end_time", "End time", Type.FLOAT),
    Definition("minimum_time_step", "Minimum time step", Type.FLOAT),
    Definition("maximum_time_step", "Maximum time step", Type.FLOAT),
    Definition("observation_time_method", "Observation time method", Type.INTEGER),
    Definition("observation_file_format", "Observation file format", Type.CHOICE, choices=("scil", "pure", "gmsh")),
    Definition("make_observation_sequence", "Make observation-time sequence", Type.BOOLEAN),
    Definition("observation_time_count", "Number of observation times", Type.INTEGER),
    Definition("observation_times", "Observation times", Type.FLOAT_LIST, count_from="observation_time_count"),
    Definition("observation_point_count", "Number of observation points", Type.INTEGER),
    Definition(
        "observation_points",
        "Observation points",
        Type.FLOAT_LIST,
        count_from="observation_point_count",
        insert_before="#define points with measured data",
    ),
    Definition("measured_point_count", "Points with measured data", Type.INTEGER),
    Definition("compute_boundary_fluxes", "Compute boundary fluxes", Type.BOOLEAN),
    Definition("print_level", "Print level", Type.INTEGER),
    Definition("nonlinear_iteration_method", "Nonlinear iteration method", Type.INTEGER),
    Definition("time_integration_method", "Time integration method", Type.INTEGER),
    Definition("inverse_modeling", "Enable inverse modeling", Type.BOOLEAN),
    Definition("integral_mass_balance", "Evaluate integral mass balance", Type.BOOLEAN),
    Definition("run_from_backup", "Run from backup", Type.BOOLEAN),
    Definition("gauss_quadrature_degree", "Gauss quadrature degree", Type.INTEGER),
)


class GlobalConfigFile(ConfigFile):
    """The DRUtES global configuration."""

    def __init__(self, path: str | Path) -> None:
        super().__init__(path, GLOBAL_DEFINITIONS)
