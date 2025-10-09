from dataclasses import dataclass
from pathlib import Path
import subprocess

import download_data

BASE = Path(__file__).parent.parent
MODELS_DIR = BASE / "data" / "models"
TARGET_EXE = BASE / "target" / "release" / "tausplit.exe"
ITERATIONS = 5
"The number of times to compute the time."


@dataclass
class ParameterSet:
    """
    A set of parameters for a test-run of the Tau-Splitting algorithm.
    """

    name: str
    network_files: list[Path]
    time: float


PARAMETER_SETS = [
    ParameterSet(
        "B-cell receptor (low)",
        [
            Path("data/models/B cell antigen receptor signaling/BCR_rxn.txt"),
            Path("data/models/B cell antigen receptor signaling/BCR_pop.txt"),
        ],
        300,
    ),
    ParameterSet(
        "B-cell receptor (high)",
        [
            Path("data/models/B cell antigen receptor signaling/BCR_rxn.txt"),
            Path("data/models/B cell antigen receptor signaling/BCR_pop_high.txt"),
        ],
        0.0009,
    ),
    ParameterSet(
        "FceRI cascade (low)",
        [
            Path("data/models/FceRI/Phosphorylation-Syk_rxn.txt"),
            Path("data/models/FceRI/Phosphorylation-Syk_pop.txt"),
        ],
        17,
    ),
    ParameterSet(
        "FceRI cascade (high)",
        [
            Path("data/models/FceRI/Phosphorylation-Syk_rxn.txt"),
            Path("data/models/FceRI/Phosphorylation-Syk_pop_high.txt"),
        ],
        0.027,
    ),
]


def main():
    download_data.main()

    subprocess.run(["cargo", "build", "--release"])

    for params in PARAMETER_SETS:
        print(f"Simulating {params.name}")
        for i in range(ITERATIONS):
            subprocess.run(
                [
                    str(x)
                    for x in [TARGET_EXE, params.time]
                    + params.network_files
                    + [
                        "--no-print-state",
                        "--count-reactions",
                        "--cpu-time",
                        "--seed",
                        str(i),
                        "--algorithm",
                        "tau-split6"
                    ]
                ],
            )


if __name__ == "__main__":
    main()
