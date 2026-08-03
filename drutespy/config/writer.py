"""Formatting-preserving writer for DRUtES configuration files."""

from __future__ import annotations

from pathlib import Path

from .parameter import Parameter, ParameterType, ScalarValue


class ConfigWriter:
    """Apply parameter values to original lines and write them atomically."""

    def render(self, lines: list[str], parameters: list[Parameter]) -> str:
        rendered = lines.copy()
        # Work bottom-up so a resized list cannot invalidate later line indexes.
        modified = [parameter for parameter in parameters if parameter.is_modified]
        ordered = sorted(
            modified,
            key=lambda item: (
                item.line_indexes[0]
                if item.line_indexes
                else item.insertion_index
                if item.insertion_index is not None
                else -1
            ),
            reverse=True,
        )
        for parameter in ordered:
            values = (
                parameter.value
                if isinstance(parameter.value, list)
                else [parameter.value]
            )
            if len(values) != len(parameter.line_indexes) and not isinstance(
                parameter.value, list
            ):
                raise ValueError(
                    f"{parameter.key!r} expects {len(parameter.line_indexes)} "
                    f"value(s), got {len(values)}."
                )
            if isinstance(parameter.value, list):
                if parameter.line_indexes:
                    first = parameter.line_indexes[0]
                    last = parameter.line_indexes[-1] + 1
                    template = rendered[first]
                elif parameter.insertion_index is not None:
                    first = last = parameter.insertion_index
                    template = "\n"
                else:
                    raise ValueError(
                        f"No insertion location is defined for {parameter.key!r}."
                    )
                rendered[first:last] = [
                    self._replace_value(
                        template, self._format(value, parameter.value_type)
                    )
                    for value in values
                ]
            else:
                index = parameter.line_indexes[0]
                rendered[index] = self._replace_value(
                    rendered[index], self._format(values[0], parameter.value_type)
                )
        return "".join(rendered)

    def write(
        self, path: Path, lines: list[str], parameters: list[Parameter]
    ) -> str:
        content = self.render(lines, parameters)
        temporary_path = path.with_name(f".{path.name}.tmp")
        temporary_path.write_text(content, encoding="utf-8", newline="")
        temporary_path.replace(path)
        return content

    @staticmethod
    def _replace_value(original: str, value: str) -> str:
        newline = "\r\n" if original.endswith("\r\n") else "\n" if original.endswith("\n") else ""
        body = original[: -len(newline)] if newline else original
        leading = body[: len(body) - len(body.lstrip())]
        trailing = body[len(body.rstrip()) :]
        return f"{leading}{value}{trailing}{newline}"

    @staticmethod
    def _format(value: ScalarValue, value_type: ParameterType) -> str:
        if value_type is ParameterType.BOOLEAN:
            return "y" if value else "n"
        return str(value)
