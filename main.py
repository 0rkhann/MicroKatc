"""Runs the hydroformylation example of the paper (ACS Catal. 2025, 15, 4739).

For your own system, write a study file and run: python microkatc.py run <study.yaml>
"""

import sys
from pathlib import Path

from microkatc import main

if __name__ == "__main__":
    example = (
        Path(__file__).resolve().parent / "examples" / "hydroformylation" / "study.yaml"
    )
    sys.exit(main(["run", str(example), *sys.argv[1:]]))
