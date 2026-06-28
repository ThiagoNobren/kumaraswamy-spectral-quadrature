from __future__ import annotations

from dataclasses import dataclass

KUMARASWAMY_CASES: list[dict[str, float | str]] = [
    {"case_id": "case_1", "label": "Kum(2,8)", "a": 2.0, "b": 8.0},
    {"case_id": "case_2", "label": "Kum(2,2.41)", "a": 2.0, "b": 2.41},
    {"case_id": "case_3", "label": "Kum(4,1.5)", "a": 4.0, "b": 1.5},
]

N_VALUES: list[int] = list(range(4, 21))

QUANTILE_LEVELS: list[float] = [0.25, 0.50, 0.75]

GRID_SIZE: int = 2000

LAMBDA_R: float = 1.0
LAMBDA_M: float = 10.0
LAMBDA_V: float = 10.0
LAMBDA_Q: float = 10.0

DOMAIN_LO: float = 0.0
DOMAIN_HI: float = 1.0


@dataclass(frozen=True)
class KumaraswamyCase:
    case_id: str
    label: str
    a: float
    b: float

    @classmethod
    def from_dict(cls, d: dict[str, float | str]) -> KumaraswamyCase:
        return cls(
            case_id=str(d["case_id"]),
            label=str(d["label"]),
            a=float(d["a"]),
            b=float(d["b"]),
        )


def all_cases() -> list[KumaraswamyCase]:
    return [KumaraswamyCase.from_dict(c) for c in KUMARASWAMY_CASES]
