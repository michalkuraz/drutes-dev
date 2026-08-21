"""Configuration parameter models."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TypeAlias

ScalarValue: TypeAlias = str | int | float | bool
ParameterValue: TypeAlias = ScalarValue | list[ScalarValue]


class ParameterType(str, Enum):
    """Types supported by the configuration editor."""

    STRING = "string"
    INTEGER = "integer"
    FLOAT = "float"
    BOOLEAN = "boolean"
    CHOICE = "choice"
    FLOAT_LIST = "float_list"
    STRING_LIST = "string_list"


@dataclass(frozen=True)
class ParameterDefinition:
    """Declarative description of a value in a positional DRUtES file."""

    key: str
    label: str
    value_type: ParameterType = ParameterType.STRING
    help_text: str = ""
    choices: tuple[str, ...] = ()
    count_from: str | None = None
    insert_before: str | None = None


@dataclass
class Parameter:
    """A parsed, editable configuration value."""

    definition: ParameterDefinition
    value: ParameterValue
    line_indexes: list[int] = field(default_factory=list)
    original_value: ParameterValue | None = None
    insertion_index: int | None = None

    def __post_init__(self) -> None:
        if self.original_value is None:
            self.original_value = (
                self.value.copy() if isinstance(self.value, list) else self.value
            )

    @property
    def is_modified(self) -> bool:
        return self.value != self.original_value

    @property
    def key(self) -> str:
        return self.definition.key

    @property
    def label(self) -> str:
        return self.definition.label

    @property
    def value_type(self) -> ParameterType:
        return self.definition.value_type

    @property
    def help_text(self) -> str:
        return self.definition.help_text

    @property
    def choices(self) -> tuple[str, ...]:
        return self.definition.choices
