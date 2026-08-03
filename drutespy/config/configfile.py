"""High-level API for loading and saving a DRUtES configuration file."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from pathlib import Path

from .parameter import Parameter, ParameterDefinition, ParameterValue
from .parser import ConfigParser
from .writer import ConfigWriter


class ConfigFile:
    """An editable configuration file backed by formatting-preserving text."""

    def __init__(
        self,
        path: str | Path,
        definitions: Sequence[ParameterDefinition],
        parser: ConfigParser | None = None,
        writer: ConfigWriter | None = None,
    ) -> None:
        self.path = Path(path)
        self.definitions = tuple(definitions)
        self.parser = parser or ConfigParser()
        self.writer = writer or ConfigWriter()
        self._lines: list[str] = []
        self._parameters: dict[str, Parameter] = {}

    def load(self) -> ConfigFile:
        """Read and parse the file, returning this instance for chaining."""
        self._lines = self.path.read_text(encoding="utf-8").splitlines(keepends=True)
        parameters = self.parser.parse(self._lines, self.definitions)
        self._parameters = {parameter.key: parameter for parameter in parameters}
        return self

    def save(self) -> None:
        """Persist current values while retaining comments and whitespace."""
        self.writer.write(self.path, self._lines, list(self._parameters.values()))
        self.load()

    def set_values(self, values: dict[str, ParameterValue]) -> None:
        """Update known parameters without coupling the model to a GUI."""
        unknown = values.keys() - self._parameters.keys()
        if unknown:
            raise KeyError(f"Unknown configuration parameter(s): {', '.join(unknown)}")
        for key, value in values.items():
            self._parameters[key].value = value

    @property
    def parameters(self) -> tuple[Parameter, ...]:
        return tuple(self._parameters.values())

    def __getitem__(self, key: str) -> Parameter:
        return self._parameters[key]

    def __iter__(self) -> Iterator[Parameter]:
        return iter(self._parameters.values())

