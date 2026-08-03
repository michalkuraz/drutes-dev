"""Parser for positional, comment-rich DRUtES configuration files."""

from __future__ import annotations

from collections.abc import Sequence

from .parameter import (
    Parameter,
    ParameterDefinition,
    ParameterType,
    ParameterValue,
)


class ConfigParseError(ValueError):
    """Raised when a configuration file does not match its definition."""


class ConfigParser:
    """Parse values while retaining their locations in the original text."""

    def parse(
        self, lines: list[str], definitions: Sequence[ParameterDefinition]
    ) -> list[Parameter]:
        value_lines = [
            index
            for index, line in enumerate(lines)
            if line.strip() and not line.lstrip().startswith("#")
        ]
        cursor = 0
        parsed: list[Parameter] = []
        by_key: dict[str, Parameter] = {}

        for definition in definitions:
            count = 1
            if definition.count_from:
                try:
                    count = int(by_key[definition.count_from].value)
                except (KeyError, TypeError, ValueError) as error:
                    raise ConfigParseError(
                        f"Invalid count parameter {definition.count_from!r} "
                        f"for {definition.key!r}."
                    ) from error

            indexes = value_lines[cursor : cursor + count]
            if len(indexes) != count:
                raise ConfigParseError(
                    f"Expected {count} value(s) for {definition.key!r}, "
                    f"but the file ended early."
                )

            raw_values = [lines[index].strip() for index in indexes]
            value: ParameterValue
            if definition.value_type is ParameterType.FLOAT_LIST:
                value = [self._convert(raw, ParameterType.FLOAT) for raw in raw_values]
            else:
                value = self._convert(raw_values[0], definition.value_type)

            insertion_index = None
            if definition.insert_before:
                insertion_index = self._find_marker(lines, definition.insert_before)
            parameter = Parameter(
                definition,
                value,
                indexes,
                insertion_index=insertion_index,
            )
            parsed.append(parameter)
            by_key[definition.key] = parameter
            cursor += count

        if cursor != len(value_lines):
            extras = ", ".join(str(index + 1) for index in value_lines[cursor:])
            raise ConfigParseError(f"Unexpected values on line(s): {extras}.")

        return parsed

    @staticmethod
    def _find_marker(lines: list[str], marker: str) -> int:
        for index, line in enumerate(lines):
            if line.strip().lower().startswith(marker.lower()):
                return index
        raise ConfigParseError(f"Required insertion marker {marker!r} was not found.")

    @staticmethod
    def _convert(raw: str, value_type: ParameterType) -> str | int | float | bool:
        try:
            if value_type is ParameterType.INTEGER:
                return int(raw)
            if value_type is ParameterType.FLOAT:
                return float(raw)
            if value_type is ParameterType.BOOLEAN:
                normalized = raw.lower()
                if normalized not in {"y", "n"}:
                    raise ValueError("expected 'y' or 'n'")
                return normalized == "y"
            return raw
        except ValueError as error:
            raise ConfigParseError(
                f"Cannot parse {raw!r} as {value_type.value}."
            ) from error
