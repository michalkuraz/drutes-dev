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

    LENGTH_UNIT_MARKER = "# GUI length unit:"

    def __init__(self, path: str | Path) -> None:
        super().__init__(path, GLOBAL_DEFINITIONS)
        self.length_unit = "m"
        self._length_unit_modified = False

    def load(self) -> GlobalConfigFile:
        super().load()
        self.length_unit = "m"
        for line in self._lines:
            if line.strip().startswith(self.LENGTH_UNIT_MARKER):
                value = line.strip()[len(self.LENGTH_UNIT_MARKER) :].strip()
                if value:
                    self.length_unit = value
                break
        self._length_unit_modified = False
        return self

    def set_length_unit(self, value: str) -> None:
        """Set GUI length units without adding a positional Fortran value."""
        self.length_unit = value
        self._length_unit_modified = True

    def save(self) -> None:
        length_unit = self.length_unit
        update_length_unit = self._length_unit_modified
        super().save()
        if not update_length_unit:
            return

        replacement = f"{self.LENGTH_UNIT_MARKER} {length_unit}\n"
        marker_index = next(
            (
                index
                for index, line in enumerate(self._lines)
                if line.strip().startswith(self.LENGTH_UNIT_MARKER)
            ),
            None,
        )
        if marker_index is None:
            if self._lines and self._lines[-1].strip():
                self._lines.append("\n")
            self._lines.append(replacement)
        else:
            self._lines[marker_index] = replacement
        temporary = self.path.with_name(f".{self.path.name}.tmp")
        temporary.write_text("".join(self._lines), encoding="utf-8", newline="")
        temporary.replace(self.path)
        self.load()
