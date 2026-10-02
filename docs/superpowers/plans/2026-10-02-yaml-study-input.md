# YAML Study Input Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `reactions.csv` + editing `main.py` with one YAML study file run by `python microkatc.py run study.yaml`, and rewrite the hydroformylation example in it with unchanged results.

**Architecture:** `study_yaml.py` (strict YAML loader) → `steps.py` (step strings) → `study.py` (`load_study` → validated `Study`) → `thermochemistry.py` (Gibbs energies and barrier tables in today's CSV format) → the existing analyses, generalized for product name, cycles, maximum yield and output step, driven by `microkatc.py` inside `<study folder>/results/`.

**Tech Stack:** Python 3.10, PyYAML 6.0.3 (new), numpy 1.26.0, pandas 1.5.3, python-copasi 4.44.295, copasi_helper `b9f26fe4`, thermochange `daed0427`, ruff 0.16.10. Tests are plain `test_*` functions in `tests/test_*.py`, each file runnable with `python tests/<file>.py` (the CI loop runs every file in its own process).

**Spec:** `docs/superpowers/specs/2026-10-02-yaml-study-input-design.md` (approved 2026-10-02, including D12–D14). Read it before any task.

## Global Constraints

- Base: `main` with #25, #26 and #27 merged (spec §0). Do not start until they are.
- Python 3.10; no new dependency besides `PyYAML==6.0.3`.
- Units live in key names (`_K`, `_M`, `_s`, `_h`, `_kcal_mol`); values are plain numbers.
- Steps: `A + B <=> C via TS`; `<=>`, `=`, `⇌` are synonyms; coefficients `2 A` or `2*A`, positive integers; every step reversible.
- Barrierless steps: TS at 4 kcal/mol above the higher side (`DIFFUSION_BARRIER = 4`, unchanged).
- File energies: thermochange `formatted_energy_outputter.sh [-g] <file> <T> <P>`, P = `0.082057366080960 * T` atm, 4th tab field × `627.509`.
- Barrier tables keep today's format and name: `G_values_of_reactions/reaction_df_{T}K_{P:.5e}atm.csv` with columns `Rx,TS,Gdir,Ginv`, `Rx` written as `A + A = B`.
- Equivalence (spec §1): barriers exact, every fitted Ea within 1e-3 kcal/mol, DRC within 1e-3.
- Exit codes: 0 success, 1 invalid study or stale results, 2 failure during the run.
- Linux/macOS only (Bash).
- Formatting: `ruff check .` and `ruff format --check .` must pass after every task.

## Review Focus

1. A studied species that is also an overall reactant (its maximum yield changes with every studied concentration) — Task 6 makes `MicroKinetics` take one threshold per concentration; tested in Task 6 and by the typed example (Task 7), whose studied species is the substrate.
2. A study run from a folder other than the repository root, or with spaces in its path — Task 7 test `test_run_from_other_folder` uses a temporary folder whose name contains a space.
3. Temperatures written as integers in YAML (`350`) versus floats in file names (`350.0`) — Task 4 converts every temperature to `float`; Task 5 test `test_barrier_table_names_use_float_temperatures` pins the file names.
4. A typed-energy study with only the microkinetics analysis (no per-temperature values needed) — Task 4 test `test_single_values_enough_without_activation_energy`.
5. Re-running after a crash mid-run (same study, partial `results/`) must reuse what exists, not refuse — Task 7 test `test_rerun_same_study_is_allowed`.

---

## File Structure

| File | Status | Responsibility |
| --- | --- | --- |
| `study_yaml.py` | create | `load_yaml(path) -> dict`: strict YAML loader (spec §3.6) |
| `steps.py` | create | `Step`, `parse_step(text)`: step and overall-reaction strings |
| `study.py` | create | `Study` dataclasses, `load_study(path)`, every check of spec §5 |
| `thermochemistry.py` | create | G per species from thermochange or typed values; `write_barrier_tables(study)` |
| `microkatc.py` | create | CLI `run`/`check`, `results/`, `run_info.json`, stale guard, runs the analyses |
| `file_operations.py` | modify | `run_directory()` replaces the import-time `CURRENT_DIRECTORY` |
| `auxiliary_functions.py` | modify | drop the G step, the cycle regex and the 1:1 yield; product-aware conversion time |
| `apparent_activation_energy.py` | modify | no G step inside; DRC takes `product` and `time_step` |
| `microkinetics_simulation.py` | modify | `MicroKinetics` takes `cycles` dict, `product`, per-concentration `max_product_M` |
| `plotting_functions.py` | modify | `PlotFunctions.grid(n)` sizes figure grids |
| `calculating_G_for_microkinetics.py` | modify | keep `GibbsEnergyCalculator` and constants; remove the command-line pipeline |
| `bash_parsing.py`, `get_G_compounds.sh`, `reactions.csv` | delete | replaced by `thermochemistry.py` and the study file |
| `GaussOutputFiles/` | move | → `examples/hydroformylation/GaussOutputFiles/` |
| `examples/hydroformylation/study.yaml` | create | the paper's example |
| `examples/typed_energies/study.yaml` | create | small cycle with typed energies |
| `main.py` | modify | one-line wrapper running the example |
| `readme_figures.py`, `tests/check_reproduction.py` | modify | read the example study and its `results/` |
| `tests/check_equivalence.py`, `tests/data/hydroformylation_baseline/` | create | compare a run with the pre-change baseline |
| `.github/workflows/nightly.yml` | modify | run the YAML example |
| `requirements.txt`, `.gitignore` | modify | PyYAML; ignore `results/` |
| `docs/USAGE.md`, `README.md`, `CLAUDE.md` | modify | document the study file |

### Task 1: Equivalence baseline from the current code

Captures today's numbers before anything changes, so Task 11 can prove the rewrite reproduces them.

**Files:**
- Create: `tests/data/hydroformylation_baseline/` (CSV and JSON, committed)
- Create: `tests/check_equivalence.py`

**Interfaces:**
- Produces: `tests/check_equivalence.py <results_dir>` — exit 0 if the run in `<results_dir>` matches the baseline within the spec §1 tolerances, exit 1 otherwise. Used by Task 11 and the nightly.

- [ ] **Step 1: Run the current pipeline in a scratch copy**

```bash
git switch main && git pull
BASE=$(git rev-parse --short HEAD)
WORK=$(mktemp -d)/run && mkdir -p "$WORK" && git archive HEAD | tar -x -C "$WORK"
(cd "$WORK" && MPLBACKEND=Agg python main.py > run.log 2>&1); echo "exit $?"
```
Expected: `exit 0` after about 3–7 minutes. `thermochange` must be exported.

- [ ] **Step 2: Save the baseline**

```bash
git switch -c yaml-study-input
OUT=tests/data/hydroformylation_baseline && mkdir -p "$OUT"
cp "$WORK"/G_values_of_reactions/reaction_df_*.csv "$OUT"/
cp "$WORK"/microkinetics_simulations/df_flux_*.csv "$OUT"/df_flux.csv
cp "$WORK"/microkinetics_simulations/df_rate_*.csv "$OUT"/df_rate.csv
cp "$WORK"/microkinetics_simulations/df_drc_*.csv "$OUT"/df_drc.csv
(cd "$WORK" && MPLBACKEND=Agg python - <<'PY') > "$OUT"/microkinetics.json
import json, numpy as np
from auxiliary_functions import AuxiliaryFunctions as A
from microkinetics_simulation import SimulationHandler
c0 = {"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 5e-4, "PMe3": 1}
conc = np.logspace(-10, -1, 19)
h = SimulationHandler(350.0, A.compute_pressure_value(350.0), "PMe3", 100_000, time_step=1)
sims = [h.get_simulation_df({**c0, "PMe3": c}) for c in conc]
cycles = A.find_intermediates_of_cycle(["0L", "1L"])
shares = A.get_concentrations_of_catalyst(sims, 1, cycles)
times = A.compute_time_of_product_conversion_given_reactant_concentration(sims, conc, "PMe3", 0.05 * 0.99)
print(json.dumps({"concentrations_M": conc.tolist(), "catalyst_M": dict(zip(cycles, shares)),
                  "t99_h": [t for t, _ in times]}, indent=1))
PY
printf 'Baseline from `main` at %s (before the YAML study change), made by Task 1 of\ndocs/superpowers/plans/2026-10-02-yaml-study-input.md. Do not regenerate.\n' "$BASE" > "$OUT"/README.md
ls "$OUT"
```
Expected: 5 `reaction_df_*.csv`, `df_flux.csv`, `df_rate.csv`, `df_drc.csv`, `microkinetics.json`, `README.md`.

- [ ] **Step 3: Write `tests/check_equivalence.py`**

```python
"""Compares a run of the example study with the baseline taken before the YAML change.

    python tests/check_equivalence.py examples/hydroformylation/results

Tolerances are the measured run-to-run noise of COPASI between two fresh runs (2026-10-02),
not targets chosen by hand:
- barrier tables: deterministic, compared exactly;
- Ea of the rate-determining steps (r08, r13): moved 1e-5 kcal/mol, limit 1e-3;
- Ea of product formation: moved 0.0044 kcal/mol, limit 0.01;
- DRC of the three key steps up to 1e-2 M PMe3: moved 5.5e-4, limit 1e-3. Above 1e-2 M less than 1 %
  of the catalyst is in the 0L cycle and the finite difference is ill-conditioned: the same input
  gave 0.9655, 0.9666 or 0.9875 for I3_1L = I4_1L at 0.1 M depending on key order, number type and
  what ran before in the process. There, and for every other step, limit 0.03 (moved up to 0.021);
- catalyst amounts at 1 h: moved 3.6e-6 relative, limit 1e-4; times to 99 %: identical, limit 1e-6.
Near-equilibrium steps (net flux a small difference of large rates) and which near-zero rows are
skipped for sign changes vary from run to run and are not compared.
"""

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "hydroformylation_baseline"
)
KEY_DRC_STEPS = ["I8_0L = I9_0L", "I3_1L = I4_1L", "I3_0L = I4_0L"]
failures = []


def same_text(series):
    return series.astype(str).str.split().str.join(" ")


def report(name, worst, limit):
    ok = worst <= limit
    print(
        f"{'ok  ' if ok else 'FAIL'} {name}: max difference {worst:.3g} (limit {limit:g})"
    )
    if not ok:
        failures.append(name)


def merged(results, kind, names):
    old = pd.read_csv(os.path.join(BASE, f"df_{kind}.csv"))
    new = pd.read_csv(
        glob.glob(
            os.path.join(results, "microkinetics_simulations", f"df_{kind}_*.csv")
        )[0]
    )
    column = next(c for c in old.columns if c.startswith("c0("))
    both = old.merge(new, on=["name", column], suffixes=("_old", "_new"))
    both = both[both["name"].isin(names)]
    expected = old[old["name"].isin(names)]
    assert len(both) == len(expected), (
        f"df_{kind}: {names} rows missing from the new run"
    )
    return both


def main(results):
    for base_file in sorted(glob.glob(os.path.join(BASE, "reaction_df_*.csv"))):
        name = os.path.basename(base_file)
        new = pd.read_csv(os.path.join(results, "G_values_of_reactions", name))
        old = pd.read_csv(base_file)
        assert (same_text(new["Rx"]) == same_text(old["Rx"])).all(), (
            f"{name}: steps differ"
        )
        difference = np.abs(
            new[["Gdir", "Ginv"]].to_numpy() - old[["Gdir", "Ginv"]].to_numpy()
        ).max()
        report(f"barriers {name}", float(difference), 1e-9)

    both = merged(results, "flux", ["r08.Flux", "r13.Flux"])
    report(
        "Ea of the rate-determining steps",
        float((both["Ea_old"] - both["Ea_new"]).abs().max()),
        1e-3,
    )
    both = merged(results, "rate", ["prod.Rate"])
    report(
        "Ea of product formation",
        float((both["Ea_old"] - both["Ea_new"]).abs().max()),
        1e-2,
    )

    def drc(path):
        return pd.read_csv(path).rename(columns=lambda c: " ".join(c.split()))

    old = drc(os.path.join(BASE, "df_drc.csv"))
    new = drc(
        glob.glob(os.path.join(results, "microkinetics_simulations", "df_drc_*.csv"))[0]
    )
    steps = [c for c in old.columns if "=" in c]
    conditioned = old[next(c for c in old.columns if c.startswith("c0("))] <= 1e-2
    key = np.abs(old[KEY_DRC_STEPS].to_numpy() - new[KEY_DRC_STEPS].to_numpy())
    report(
        "DRC of the key steps up to 1e-2 M",
        float(key[conditioned.to_numpy()].max()),
        1e-3,
    )
    report(
        "DRC of the key steps above 1e-2 M",
        float(key[~conditioned.to_numpy()].max(initial=0)),
        3e-2,
    )
    report(
        "DRC of every step",
        float(np.abs(old[steps].to_numpy() - new[steps].to_numpy()).max()),
        3e-2,
    )

    with open(os.path.join(BASE, "microkinetics.json")) as f:
        old = json.load(f)
    with open(os.path.join(results, "microkinetics.json")) as f:
        new = json.load(f)
    for cycle in old["catalyst_M"]:
        a, b = np.array(old["catalyst_M"][cycle]), np.array(new["catalyst_M"][cycle])
        report(
            f"catalyst in {cycle} at 1 h",
            float(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30))),
            1e-4,
        )
    a, b = np.array(old["t99_h"]), np.array(new["t99_h"])
    report("time to 99 % conversion", float(np.max(np.abs(a - b) / a)), 1e-6)

    if failures:
        sys.exit(f"{len(failures)} comparison(s) failed: {', '.join(failures)}")
    print("Run matches the baseline.")


if __name__ == "__main__":
    main(sys.argv[1])
```

`microkinetics.json` in the results folder is written by `microkatc.py` (Task 7) in the same format.

- [ ] **Step 4: Commit**

```bash
ruff check tests/check_equivalence.py && ruff format tests/check_equivalence.py
git add tests/data/hydroformylation_baseline tests/check_equivalence.py
git commit -m "test: baseline of the hydroformylation results before the YAML study change"
```

---

### Task 2: Strict YAML loader

**Files:**
- Create: `study_yaml.py`
- Create: `tests/test_study_yaml.py`
- Modify: `requirements.txt` (add `PyYAML==6.0.3`)

**Interfaces:**
- Produces: `load_yaml(path: str | os.PathLike) -> dict`; raises `StudyError(messages: list[str])`.
- Produces: `class StudyError(Exception)` with attribute `messages: list[str]`; `str(error)` joins them with newlines. Defined here, imported by `steps.py`, `study.py`, `microkatc.py`.

- [ ] **Step 1: Write the failing tests** (`tests/test_study_yaml.py`)

```python
"""Strict YAML loading for study files (spec section 3.6).

python tests/test_study_yaml.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from study_yaml import StudyError, load_yaml


def write(text):
    path = os.path.join(tempfile.mkdtemp(), "study.yaml")
    with open(path, "w") as f:
        f.write(text)
    return path


def test_exponent_numbers_without_a_dot_are_numbers():
    doc = load_yaml(write("a: 1e-10\nb: 5e-4\nc: 1.0e-10\nd: 19\ne: 0.1\n"))
    assert doc == {"a": 1e-10, "b": 5e-4, "c": 1e-10, "d": 19, "e": 0.1}
    assert isinstance(doc["d"], int)


def test_species_names_like_no_stay_text():
    doc = load_yaml(
        write("initial_M: {NO: 0.05, ON: 1, Y: 2, yes: 3, off: 4}\nflag: true\n")
    )
    assert doc["initial_M"] == {"NO": 0.05, "ON": 1, "Y": 2, "yes": 3, "off": 4}
    assert doc["flag"] is True


def test_temperature_keys_are_numbers():
    doc = load_yaml(write("TS1: {325: 11.5, 337.5: 11.6}\n"))
    assert doc["TS1"] == {325: 11.5, 337.5: 11.6}


def test_duplicate_key_is_an_error_with_both_lines():
    try:
        load_yaml(
            write("species:\n  energies_kcal_mol:\n    TS1: 1\n    A: 2\n    TS1: 3\n")
        )
    except StudyError as error:
        assert "TS1 appears twice (lines 3 and 5)" in str(error), str(error)
    else:
        raise AssertionError("expected StudyError")


def test_syntax_error_names_the_line():
    try:
        load_yaml(write("a: 1\nb: [1, 2\nc: 3\n"))
    except StudyError as error:
        assert "line" in str(error) and "study.yaml" in str(error), str(error)
    else:
        raise AssertionError("expected StudyError")


def test_top_level_must_be_a_mapping():
    try:
        load_yaml(write("- a\n- b\n"))
    except StudyError as error:
        assert "must be a mapping" in str(error)
    else:
        raise AssertionError("expected StudyError")


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
```

- [ ] **Step 2: Run to see it fail**

Run: `python tests/test_study_yaml.py`
Expected: `ModuleNotFoundError: No module named 'study_yaml'`

- [ ] **Step 3: Implement `study_yaml.py`**

```python
"""Strict YAML loading for MicroKatc study files.

PyYAML follows YAML 1.1, which reads 1e-10 as text, NO/ON/yes/off as booleans, and keeps the
last of two duplicate keys silently. Study files need the opposite on all three counts.
"""

import os
import re

import yaml


class StudyError(Exception):
    """One or more problems in a study file; messages lists each with its location"""

    def __init__(self, messages):
        self.messages = list(messages)
        super().__init__("\n".join(self.messages))


class StrictLoader(yaml.SafeLoader):
    """SafeLoader with YAML 1.2-style floats, lowercase-only booleans and duplicate-key errors"""

    def construct_mapping(self, node, deep=False):
        lines = {}
        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)
            line = key_node.start_mark.line + 1
            if key in lines:
                raise yaml.constructor.ConstructorError(
                    None,
                    None,
                    f"{key} appears twice (lines {lines[key]} and {line})",
                    key_node.start_mark,
                )
            lines[key] = line
        return super().construct_mapping(node, deep=deep)


# Copy the resolver table so SafeLoader itself is not changed
StrictLoader.yaml_implicit_resolvers = {
    first: [
        (tag, regexp)
        for tag, regexp in resolvers
        if tag not in ("tag:yaml.org,2002:bool", "tag:yaml.org,2002:float")
    ]
    for first, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
StrictLoader.add_implicit_resolver(
    "tag:yaml.org,2002:bool", re.compile(r"^(?:true|false)$"), list("tf")
)
StrictLoader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(
        r"""^[-+]?(?:(?:[0-9][0-9_]*\.[0-9_]*|\.[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)?
        |[0-9][0-9_]*[eE][-+]?[0-9]+
        |\.(?:inf|Inf|INF)
        |\.(?:nan|NaN|NAN))$""",
        re.VERBOSE,
    ),
    list("-+0123456789."),
)


def load_yaml(path):
    """Reads a study file; raises StudyError with the file name and line for any YAML problem"""
    name = os.fspath(path)
    try:
        with open(name, encoding="utf-8") as f:
            doc = yaml.load(f, Loader=StrictLoader)
    except yaml.YAMLError as error:
        mark = getattr(error, "problem_mark", None) or getattr(
            error, "context_mark", None
        )
        where = f"{name} line {mark.line + 1}" if mark else name
        problem = getattr(error, "problem", None) or str(error)
        raise StudyError([f"{where}: {problem}"]) from error
    if not isinstance(doc, dict):
        raise StudyError(
            [
                f"{name}: the study must be a mapping of keys such as species, steps and conditions"
            ]
        )
    return doc
```

- [ ] **Step 4: Add PyYAML and run the tests**

Add the line `PyYAML==6.0.3` to `requirements.txt` after `pandas==1.5.3`, then:

```bash
pip install PyYAML==6.0.3
python tests/test_study_yaml.py
```
Expected: `ok`

- [ ] **Step 5: Commit**

```bash
ruff check . && ruff format --check .
git add study_yaml.py tests/test_study_yaml.py requirements.txt
git commit -m "feat: strict YAML loader for study files"
```

---

### Task 3: Step strings

**Files:**
- Create: `steps.py`
- Create: `tests/test_steps.py`

**Interfaces:**
- Consumes: `StudyError` is not used here; `parse_step` raises `StepSyntaxError(ValueError)`, which `study.py` turns into located messages.
- Produces:
  - `@dataclass(frozen=True) class Step: reactants: tuple[tuple[str, int], ...]; products: tuple[tuple[str, int], ...]; ts: str | None; text: str`
  - `Step.species() -> list[str]` — reactants then products, first-seen order, no TS
  - `Step.copasi_equation() -> str` — coefficients expanded, e.g. `"A + A = B"`
  - `Step.coefficient(name: str, side: str) -> int` — `side` is `"reactants"` or `"products"`; 0 if absent
  - `parse_step(text: str, allow_ts: bool = True) -> Step`
  - `class StepSyntaxError(ValueError)`

- [ ] **Step 1: Write the failing tests** (`tests/test_steps.py`)

```python
"""Parsing of step and overall-reaction strings (spec section 3.3).

python tests/test_steps.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from steps import StepSyntaxError, parse_step


def test_simple_barrierless_step():
    step = parse_step("I1_0L <=> I2_0L + CO")
    assert step.reactants == (("I1_0L", 1),)
    assert step.products == (("I2_0L", 1), ("CO", 1))
    assert step.ts is None
    assert step.copasi_equation() == "I1_0L = I2_0L + CO"


def test_step_with_transition_state():
    step = parse_step("I6_0L + H2 <=> I8_0L  via TS3_0L")
    assert step.ts == "TS3_0L"
    assert step.species() == ["I6_0L", "H2", "I8_0L"]


def test_separators_are_synonyms():
    for text in ("A + B <=> C", "A + B = C", "A + B ⇌ C", "A+B<=>C"):
        assert parse_step(text).copasi_equation() == "A + B = C", text


def test_coefficients_are_expanded():
    for text in ("2 A <=> A2 via TS2", "2*A <=> A2 via TS2", "A + A <=> A2 via TS2"):
        step = parse_step(text)
        assert step.coefficient("A", "reactants") == 2, text
        assert step.copasi_equation() == "A + A = A2", text


def test_errors():
    cases = {
        "I4_0L + CO I5_0L": "no <=> between the two sides",
        "A <=> B <=> C": "more than one <=>",
        "A <=> ": "the right side is empty",
        "1.5 A <=> B": "coefficients must be positive integers",
        "0 A <=> B": "coefficients must be positive integers",
        "A + <=> B": "empty species name",
        "A <=> B via": "via must be followed by one transition-state name",
        "A <=> B via TS1 TS2": "via must be followed by one transition-state name",
        "A b <=> B": "names cannot contain spaces",
        "1A <=> B": "names must start with a letter",
    }
    for text, message in cases.items():
        try:
            parse_step(text)
        except StepSyntaxError as error:
            assert message in str(error), (text, str(error))
        else:
            raise AssertionError(f"no error for {text!r}")


def test_overall_reaction_rejects_via():
    try:
        parse_step("ete + CO + H2 <=> prod via TS1", allow_ts=False)
    except StepSyntaxError as error:
        assert "via is not allowed here" in str(error)
    else:
        raise AssertionError("expected StepSyntaxError")


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
```

- [ ] **Step 2: Run to see it fail**

Run: `python tests/test_steps.py`
Expected: `ModuleNotFoundError: No module named 'steps'`

- [ ] **Step 3: Implement `steps.py`**

```python
"""Step strings of a study: "A + B <=> C via TS", with optional integer coefficients"""

import re
from dataclasses import dataclass

SEPARATOR = re.compile(r"<=>|⇌|=")
NAME = re.compile(r"[A-Za-z][^\s+*=⇌<>]*")
TERM = re.compile(r"^(?:(?P<coef>\d+)\s*\*?\s+|(?P<coef2>\d+)\*)?(?P<name>.*)$")


class StepSyntaxError(ValueError):
    """A step or overall-reaction string that cannot be parsed"""


@dataclass(frozen=True)
class Step:
    """One reversible elementary step; coefficients are kept per species"""

    reactants: tuple
    products: tuple
    ts: object
    text: str

    def species(self):
        """Reactant and product names, first-seen order, without the transition state"""
        seen = []
        for name, _ in self.reactants + self.products:
            if name not in seen:
                seen.append(name)
        return seen

    def coefficient(self, name, side):
        """Coefficient of name on side ("reactants" or "products"), 0 if absent"""
        return dict(getattr(self, side)).get(name, 0)

    def copasi_equation(self):
        """The step as copasi_helper reads it, coefficients written as repeated species"""

        def side(terms):
            return " + ".join(name for name, n in terms for _ in range(n))

        return f"{side(self.reactants)} = {side(self.products)}"


def _side(text, label, whole):
    if not text.strip():
        raise StepSyntaxError(f'"{whole}": the {label} side is empty')
    counts = {}
    for raw in text.split("+"):
        term = raw.strip()
        if not term:
            raise StepSyntaxError(f'"{whole}": empty species name next to a +')
        match = TERM.match(term)
        coef = match.group("coef") or match.group("coef2")
        name = match.group("name").strip()
        if (
            re.match(r"^\d*\.\d+|^0+\s*\*?\s", term)
            or coef == "0"
            or (coef is not None and int(coef) == 0)
        ):
            raise StepSyntaxError(
                f'"{whole}": coefficients must be positive integers ({term})'
            )
        if " " in name:
            raise StepSyntaxError(f'"{whole}": names cannot contain spaces ({name})')
        if not NAME.fullmatch(name):
            raise StepSyntaxError(f'"{whole}": names must start with a letter ({name})')
        counts[name] = counts.get(name, 0) + (int(coef) if coef else 1)
    return tuple(counts.items())


def parse_step(text, allow_ts=True):
    """Parses "A + 2 B <=> C via TS"; raises StepSyntaxError with the reason"""
    whole = " ".join(str(text).split())
    body, ts = whole, None
    if re.search(r"\svia(\s|$)", whole):
        if not allow_ts:
            raise StepSyntaxError(f'"{whole}": via is not allowed here')
        body, _, after = whole.partition(" via")
        names = after.split()
        if len(names) != 1:
            raise StepSyntaxError(
                f'"{whole}": via must be followed by one transition-state name'
            )
        ts = names[0]
    sides = SEPARATOR.split(body)
    if len(sides) == 1:
        raise StepSyntaxError(f'"{whole}": no <=> between the two sides')
    if len(sides) > 2:
        raise StepSyntaxError(f'"{whole}": more than one <=>')
    left, right = sides
    return Step(_side(left, "left", whole), _side(right, "right", whole), ts, whole)
```

- [ ] **Step 4: Run the tests**

Run: `python tests/test_steps.py`
Expected: `ok`. If a case in `test_errors` fails, fix `_side`, not the test: each message is promised by spec §5.

- [ ] **Step 5: Commit**

```bash
ruff check . && ruff format --check .
git add steps.py tests/test_steps.py
git commit -m "feat: parse study step strings with coefficients and inline transition states"
```

---

### Task 4: `load_study` and every input check

**Files:**
- Create: `study.py`
- Create: `tests/test_study.py`

**Interfaces:**
- Consumes: `load_yaml`, `StudyError` (Task 2); `parse_step`, `StepSyntaxError`, `Step` (Task 3).
- Produces:
  - `load_study(path, check_energy_sources=True) -> Study`; raises `StudyError` with every problem. `check_energy_sources=False` skips output-file and `$thermochange` checks (for scripts reading finished results).
  - `@dataclass(frozen=True) class Study` with fields `path: Path`, `name: str`, `files: Path | None`, `energies_kcal_mol: dict[str, dict[float, float]] | None`, `product: str`, `overall: Step`, `cycles: dict[str, tuple[str, ...]]`, `vibrational_correction: str`, `steps: tuple[Step, ...]`, `temperature_K: float`, `studied_species: str`, `studied_concentrations_M: tuple[float, ...]`, `output_step_s: float`, `activation_energy: ActivationEnergy | None`, `degree_of_rate_control: DegreeOfRateControl | None`, `microkinetics: Microkinetics | None`, `warnings: tuple[str, ...]`.
  - Methods: `results_dir -> Path` (property, `<study folder>/results`), `species() -> list[str]`, `transition_states() -> list[str]`, `temperatures() -> list[float]`, `reactions() -> list[str]` (copasi equations in step order), `max_product_M(initial_M: dict) -> float`.
  - `ActivationEnergy(initial_M, temperatures_K: tuple, simulation_time_s, sampling_time_h)`, `DegreeOfRateControl(initial_M, e_shift_kcal_mol, simulation_time_s, cores: int)`, `Microkinetics(initial_M, simulation_time_s, catalyst_snapshot_h, conversion, plot_species: tuple)`.
  - Every temperature is a `float` (Review Focus 3); a single typed energy is stored as `{temperature_K: value}`.

- [ ] **Step 1: Write the failing tests** (`tests/test_study.py`)

```python
"""Loading and checking study files (spec sections 3 and 5).

python tests/test_study.py
"""

import copy
import os
import sys
import tempfile

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from study import load_study
from study_yaml import StudyError

TYPED = {
    "name": "Typed toy",
    "species": {
        "energies_kcal_mol": {"S": 0.0, "P": -10.0, "C1": 0.0, "C2": -2.0, "TS1": 15.0},
        "product": "P",
        "overall_reaction": "S <=> P",
        "cycles": {"cat": ["C1", "C2"]},
    },
    "steps": ["C1 + S <=> C2", "C2 <=> C1 + P  via TS1"],
    "conditions": {
        "temperature_K": 300,
        "studied_species": "S",
        "studied_range_M": {"from": 0.01, "to": 0.1, "points": 3},
    },
    "analyses": {
        "microkinetics": {
            "initial_M": {"C1": 1e-3},
            "simulation_time_s": 1000,
            "catalyst_snapshot_h": 0.01,
            "conversion": 0.99,
            "plot_species": ["C1", "C2"],
        }
    },
}


def write(doc, extra_files=()):
    folder = tempfile.mkdtemp()
    for name in extra_files:
        os.makedirs(os.path.join(folder, os.path.dirname(name)), exist_ok=True)
        open(os.path.join(folder, name), "w").close()
    path = os.path.join(folder, "study.yaml")
    with open(path, "w") as f:
        yaml.safe_dump(doc, f, allow_unicode=True)
    return path


def errors_of(doc, **kwargs):
    try:
        load_study(write(doc, **kwargs))
    except StudyError as error:
        return error.messages
    raise AssertionError("expected StudyError")


def changed(edit):
    doc = copy.deepcopy(TYPED)
    edit(doc)
    return doc


def assert_error(doc, text, **kwargs):
    messages = errors_of(doc, **kwargs)
    assert any(text in m for m in messages), (text, messages)


def test_typed_study_loads_with_derived_values():
    study = load_study(write(TYPED))
    assert study.temperature_K == 300.0 and isinstance(study.temperature_K, float)
    assert (
        len(study.studied_concentrations_M) == 3
        and study.studied_concentrations_M[0] == 0.01
    )
    assert abs(study.studied_concentrations_M[1] - 10**-1.5) < 1e-12
    assert study.species() == ["C1", "S", "C2", "P"]
    assert study.transition_states() == ["TS1"]
    assert study.energies_kcal_mol["TS1"] == {300.0: 15.0}
    assert study.reactions() == ["C1 + S = C2", "C2 = C1 + P"]
    assert study.cycles == {"cat": ("C1", "C2")}
    assert study.results_dir.name == "results"
    assert study.degree_of_rate_control is None and study.activation_energy is None


def test_max_product_uses_coefficients():
    doc = changed(lambda d: d["species"].update(overall_reaction="2 S + C1 <=> P"))
    doc["species"]["cycles"] = {"cat": ["C2"]}
    doc["analyses"]["microkinetics"]["initial_M"] = {"C2": 1e-3}
    study = load_study(write(doc))
    assert study.max_product_M({"S": 1.0, "C1": 0.3}) == 0.3
    assert study.max_product_M({"S": 0.4, "C1": 0.3}) == 0.2


def test_unknown_key_suggests_the_closest():
    assert_error(
        changed(lambda d: d["conditions"].update(temprature_K=300)),
        "conditions.temprature_K: unknown key; did you mean temperature_K?",
    )


def test_files_and_energies_together():
    assert_error(
        changed(lambda d: d["species"].update(files="out")),
        "species: give files or energies_kcal_mol, not both",
    )


def test_required_keys():
    assert_error(
        changed(lambda d: d["species"].pop("product")), "species.product: required"
    )
    assert_error(
        changed(lambda d: d["species"].pop("overall_reaction")),
        "species.overall_reaction: required",
    )
    assert_error(changed(lambda d: d.pop("analyses")), "analyses: required")
    assert_error(
        changed(lambda d: d.update(analyses={})), "analyses: request at least one"
    )


def test_step_errors_are_located():
    assert_error(
        changed(lambda d: d["steps"].append("C2 + S C1")),
        'steps[2] "C2 + S C1": no <=> between the two sides',
    )
    assert_error(changed(lambda d: d["steps"].append("1.5 S <=> P")), "steps[2]")
    assert_error(
        changed(lambda d: d["steps"].append("C2 <=> P + C1 via TS1")),
        "steps[2] repeats steps[1]",
    )


def test_missing_energy_source():
    assert_error(
        changed(lambda d: d["species"]["energies_kcal_mol"].pop("TS1")),
        "no energy for TS1: add it to species.energies_kcal_mol",
    )
    doc = changed(lambda d: d["species"].pop("energies_kcal_mol"))
    doc["species"]["files"] = "out"
    os.environ["thermochange"] = tempfile.mkdtemp()
    messages = errors_of(
        doc, extra_files=["out/S.out", "out/P.out", "out/C1.out", "out/C2.out"]
    )
    assert "no energy for TS1: out/TS1.out not found" in messages, messages
    assert any(
        "does not contain formatters/formatted_energy_outputter.sh" in m
        for m in messages
    ), messages
    del os.environ["thermochange"]
    assert_error(doc, "$thermochange is not set", extra_files=["out/TS1.out"])
    # Scripts that only read finished results skip these checks
    from study import load_study as load

    study = load(write(doc), check_energy_sources=False)
    assert study.files.name == "out"


def test_name_consistency():
    assert_error(
        changed(lambda d: d["species"]["cycles"]["cat"].append("C9")),
        "species.cycles.cat: C9 does not appear in any step",
    )
    assert_error(
        changed(lambda d: d["species"]["cycles"].update(other=["C2"])),
        "species.cycles.other: C2 is also in cycle cat",
    )
    assert_error(
        changed(lambda d: d["species"].update(product="Q")),
        "species.product: Q does not appear in any step",
    )
    assert_error(
        changed(lambda d: d["species"].update(overall_reaction="S + X <=> P")),
        "species.overall_reaction: X does not appear in any step",
    )
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"].update(plot_species=["Z"])),
        "plot_species: Z does not appear in any step",
    )


def test_studied_species_and_product_not_in_initial_M():
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"]["initial_M"].update(S=0.1)),
        "the studied species is set by conditions.studied_range_M",
    )
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"]["initial_M"].update(P=0.0)),
        "the product always starts at 0",
    )
    assert_error(
        changed(
            lambda d: d["analyses"]["microkinetics"].update(
                initial_M={"S": 0.0, "P": 0.0}
            )
        ),
        "at least one catalyst intermediate",
    )


def test_typed_energies_need_every_activation_energy_temperature():
    def add_ea(d):
        d["analyses"]["activation_energy"] = {
            "initial_M": {"C1": 1e-6},
            "temperatures_K": {"from": 290, "to": 310, "points": 3},
            "simulation_time_s": 1000,
            "sampling_time_h": 0.1,
        }

    assert_error(
        changed(add_ea),
        "species.energies_kcal_mol.TS1: needs values at 290, 310 K for activation_energy",
    )


def test_single_values_enough_without_activation_energy():
    study = load_study(write(TYPED))
    assert study.temperatures() == [300.0]


def test_numbers():
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"].update(conversion=1.2)),
        "conversion: must be between 0 and 1",
    )
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"].update(catalyst_snapshot_h=1)),
        "1 h is after the end of the simulation (1000 s = 0.28 h)",
    )
    assert_error(
        changed(lambda d: d["conditions"].update(output_step_s=7)),
        "0.01 h is not a multiple of output_step_s (7 s)",
    )
    assert_error(
        changed(lambda d: d["conditions"].update(temperature_K=-5)),
        "conditions.temperature_K: must be positive",
    )
    assert_error(
        changed(
            lambda d: d["conditions"]["studied_range_M"].update({"from": 1, "to": 0.1})
        ),
        "from must be smaller than to",
    )


def test_unassigned_species_warning():
    doc = changed(lambda d: d["steps"].append("C2 + L <=> C3"))
    doc["species"]["energies_kcal_mol"].update(L=0.0, C3=-1.0)
    study = load_study(write(doc))
    assert study.warnings == (
        "warning: L, C3 belong to no cycle, so the catalyst distribution does not count them",
    ), study.warnings


def test_all_errors_reported_together():
    doc = changed(lambda d: d["conditions"].update(temprature_K=1))
    doc["species"]["product"] = "Q"
    messages = errors_of(doc)
    assert len(messages) >= 2, messages


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
```

- [ ] **Step 2: Run to see it fail**

Run: `python tests/test_study.py`
Expected: `ModuleNotFoundError: No module named 'study'`

- [ ] **Step 3: Implement `study.py`**

```python
"""A MicroKatc study: the YAML file read, checked and turned into one Study object (spec sections 3 and 5)"""

import difflib
import math
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from steps import StepSyntaxError, parse_step
from study_yaml import StudyError, load_yaml

RANGE = {"from": None, "to": None, "points": None}
SCHEMA = {
    "name": None,
    "species": {
        "files": None,
        "energies_kcal_mol": None,
        "product": None,
        "overall_reaction": None,
        "cycles": None,
        "vibrational_correction": None,
    },
    "steps": None,
    "conditions": {
        "temperature_K": None,
        "studied_species": None,
        "studied_range_M": RANGE,
        "output_step_s": None,
    },
    "analyses": {
        "activation_energy": {
            "initial_M": None,
            "temperatures_K": RANGE,
            "simulation_time_s": None,
            "sampling_time_h": None,
        },
        "degree_of_rate_control": {
            "initial_M": None,
            "e_shift_kcal_mol": None,
            "simulation_time_s": None,
            "cores": None,
        },
        "microkinetics": {
            "initial_M": None,
            "simulation_time_s": None,
            "catalyst_snapshot_h": None,
            "conversion": None,
            "plot_species": None,
        },
    },
}
THERMOCHANGE_SCRIPT = os.path.join("formatters", "formatted_energy_outputter.sh")


@dataclass(frozen=True)
class ActivationEnergy:
    initial_M: dict
    temperatures_K: tuple
    simulation_time_s: float
    sampling_time_h: float


@dataclass(frozen=True)
class DegreeOfRateControl:
    initial_M: dict
    e_shift_kcal_mol: float
    simulation_time_s: float
    cores: int


@dataclass(frozen=True)
class Microkinetics:
    initial_M: dict
    simulation_time_s: float
    catalyst_snapshot_h: float
    conversion: float
    plot_species: tuple


@dataclass(frozen=True)
class Study:
    """Everything the analyses need, validated; built only by load_study"""

    path: Path
    name: str
    files: object  # Path to the output-file folder, or None with typed energies
    energies_kcal_mol: object  # {species: {temperature_K: G}}, or None with files
    product: str
    overall: object  # steps.Step
    cycles: dict
    vibrational_correction: str
    steps: tuple
    temperature_K: float
    studied_species: str
    studied_concentrations_M: tuple
    output_step_s: float
    activation_energy: object = None
    degree_of_rate_control: object = None
    microkinetics: object = None
    warnings: tuple = field(default_factory=tuple)

    @property
    def results_dir(self):
        return self.path.parent / "results"

    def species(self):
        """Every species in the steps, first-seen order (transition states excluded)"""
        seen = []
        for step in self.steps:
            seen += [name for name in step.species() if name not in seen]
        return seen

    def transition_states(self):
        return [step.ts for step in self.steps if step.ts is not None]

    def temperatures(self):
        """Every temperature at which Gibbs energies are needed, as floats"""
        needed = {self.temperature_K}
        if self.activation_energy is not None:
            needed.update(self.activation_energy.temperatures_K)
        return sorted(needed)

    def reactions(self):
        """The steps as copasi_helper equations, in file order (COPASI names them r01, r02, ...)"""
        return [step.copasi_equation() for step in self.steps]

    def max_product_M(self, initial_M):
        """Largest product concentration the overall reaction allows from these initial concentrations"""
        limiting = min(
            initial_M.get(name, 0.0) / n for name, n in self.overall.reactants
        )
        return limiting * self.overall.coefficient(self.product, "products")


def _unknown_keys(doc, schema, prefix, errors):
    for key, value in doc.items():
        where = f"{prefix}{key}"
        if key not in schema:
            close = difflib.get_close_matches(str(key), [str(k) for k in schema], n=1)
            hint = f"; did you mean {close[0]}?" if close else ""
            errors.append(f"{where}: unknown key{hint}")
        elif isinstance(schema[key], dict) and isinstance(value, dict):
            _unknown_keys(value, schema[key], f"{where}.", errors)


def _number(value, where, errors, positive=True, integer=False):
    ok_type = isinstance(value, int) if integer else isinstance(value, (int, float))
    if (
        isinstance(value, bool)
        or not ok_type
        or (isinstance(value, float) and not math.isfinite(value))
    ):
        errors.append(
            f"{where}: must be {'an integer' if integer else 'a number'}, got {value!r}"
        )
        return None
    if positive and value <= 0:
        errors.append(f"{where}: must be positive, got {value}")
        return None
    return float(value) if not integer else int(value)


def _range(value, where, errors, minimum_points=1, log=False):
    if not isinstance(value, dict) or set(value) != {"from", "to", "points"}:
        errors.append(f"{where}: must be {{from: ..., to: ..., points: ...}}")
        return None
    start = _number(value["from"], f"{where}.from", errors)
    stop = _number(value["to"], f"{where}.to", errors)
    points = _number(value["points"], f"{where}.points", errors, integer=True)
    if None in (start, stop, points):
        return None
    if points < minimum_points:
        errors.append(f"{where}.points: needs at least {minimum_points}, got {points}")
        return None
    if points > 1 and not start < stop:
        errors.append(f"{where}: from must be smaller than to")
        return None
    if log:
        return tuple(
            float(v) for v in np.logspace(np.log10(start), np.log10(stop), points)
        )
    return tuple(float(v) for v in np.linspace(start, stop, points))


def _concentrations(value, where, errors):
    if not isinstance(value, dict) or not value:
        errors.append(f"{where}: must map species to initial concentrations in M")
        return {}
    out = {}
    for name, c in value.items():
        number = _number(c, f"{where}.{name}", errors, positive=False)
        if number is not None and number < 0:
            errors.append(f"{where}.{name}: must not be negative")
        elif number is not None:
            out[str(name)] = number
    return out


def _on_grid(hours, step_s, total_s, where, errors):
    seconds = hours * 3600
    if seconds > total_s:
        errors.append(
            f"{where}: {hours:g} h is after the end of the simulation ({total_s:g} s = {total_s / 3600:.2g} h)"
        )
    elif abs(seconds / step_s - round(seconds / step_s)) > 1e-9:
        errors.append(
            f"{where}: {hours:g} h is not a multiple of output_step_s ({step_s:g} s)"
        )


def _energies(
    value, where, errors, needed_temperatures, working_temperature, ea_requested
):
    """{species: {T: G}}; a single number is valid at the working temperature only"""
    if not isinstance(value, dict) or not value:
        errors.append(f"{where}: must map species to Gibbs energies in kcal/mol")
        return {}
    out = {}
    for name, entry in value.items():
        if isinstance(entry, dict):
            table = {}
            for temperature, g in entry.items():
                t = _number(temperature, f"{where}.{name} temperature", errors)
                v = _number(g, f"{where}.{name}.{temperature}", errors, positive=False)
                if t is not None and v is not None:
                    table[t] = v
        else:
            v = _number(entry, f"{where}.{name}", errors, positive=False)
            table = {working_temperature: v} if v is not None else {}
        missing = [
            t
            for t in needed_temperatures
            if not any(abs(t - known) < 1e-6 for known in table)
        ]
        if missing and table:
            purpose = "activation_energy" if ea_requested else "the working temperature"
            listed = ", ".join(f"{t:g}" for t in missing)
            errors.append(f"{where}.{name}: needs values at {listed} K for {purpose}")
        out[str(name)] = table
    return out


def load_study(path, check_energy_sources=True):
    """Reads and checks a study file; raises StudyError listing every problem found.

    check_energy_sources=False skips the output-file and $thermochange checks, for scripts that
    only read the results of a finished run.
    """
    path = Path(path).resolve()
    doc = load_yaml(path)
    errors, warnings = [], []
    _unknown_keys(doc, SCHEMA, "", errors)

    species_block = doc.get("species") if isinstance(doc.get("species"), dict) else {}
    conditions = (
        doc.get("conditions") if isinstance(doc.get("conditions"), dict) else {}
    )
    analyses = doc.get("analyses") if isinstance(doc.get("analyses"), dict) else {}
    for key in ("species", "steps", "conditions", "analyses"):
        if key not in doc:
            errors.append(f"{key}: required")
    for key in ("product", "overall_reaction"):
        if key not in species_block:
            errors.append(f"species.{key}: required")
    if ("files" in species_block) == ("energies_kcal_mol" in species_block):
        errors.append(
            "species: give files or energies_kcal_mol, not both"
            if "files" in species_block
            else "species: give files or energies_kcal_mol"
        )
    if analyses is not None and not any(k in analyses for k in SCHEMA["analyses"]):
        errors.append(
            "analyses: request at least one of activation_energy, degree_of_rate_control, microkinetics"
        )

    # Steps
    steps, raw_steps = [], doc.get("steps", [])
    if not isinstance(raw_steps, list) or not raw_steps:
        errors.append("steps: must be a list of steps such as 'A + B <=> C via TS'")
        raw_steps = []
    seen = {}
    for i, text in enumerate(raw_steps):
        try:
            step = parse_step(text)
        except StepSyntaxError as error:
            errors.append(f"steps[{i}] {error}")
            continue
        key = (tuple(sorted(step.reactants)), tuple(sorted(step.products)), step.ts)
        if key in seen:
            errors.append(f"steps[{i}] repeats steps[{seen[key]}] ({step.text})")
        seen.setdefault(key, i)
        steps.append(step)
    names = []
    for step in steps:
        names += [n for n in step.species() if n not in names]
    transition_states = [s.ts for s in steps if s.ts is not None]

    # Conditions
    temperature = (
        _number(conditions.get("temperature_K"), "conditions.temperature_K", errors)
        if "temperature_K" in conditions
        else None
    )
    if "temperature_K" not in conditions:
        errors.append("conditions.temperature_K: required")
    studied = conditions.get("studied_species")
    if studied is None:
        errors.append("conditions.studied_species: required")
    elif names and studied not in names:
        errors.append(
            f"conditions.studied_species: {studied} does not appear in any step"
        )
    concentrations = None
    if "studied_range_M" in conditions:
        concentrations = _range(
            conditions["studied_range_M"],
            "conditions.studied_range_M",
            errors,
            log=True,
        )
    else:
        errors.append("conditions.studied_range_M: required")
    output_step = _number(
        conditions.get("output_step_s", 1), "conditions.output_step_s", errors
    )

    # Species block
    product = species_block.get("product")
    if product is not None and names and product not in names:
        errors.append(f"species.product: {product} does not appear in any step")
    overall = None
    if "overall_reaction" in species_block:
        try:
            overall = parse_step(species_block["overall_reaction"], allow_ts=False)
        except StepSyntaxError as error:
            errors.append(f"species.overall_reaction {error}")
    if overall is not None:
        for name in overall.species():
            if names and name not in names:
                errors.append(
                    f"species.overall_reaction: {name} does not appear in any step"
                )
        if product is not None and overall.coefficient(product, "products") == 0:
            errors.append(
                f"species.overall_reaction: the product {product} must be on its right side"
            )
    cycles, raw_cycles = {}, species_block.get("cycles")
    if raw_cycles is not None:
        if not isinstance(raw_cycles, dict) or not all(
            isinstance(v, list) for v in raw_cycles.values()
        ):
            errors.append(
                "species.cycles: must map each cycle label to a list of intermediates"
            )
        else:
            owner = {}
            for label, members in raw_cycles.items():
                for member in members:
                    if names and member not in names:
                        errors.append(
                            f"species.cycles.{label}: {member} does not appear in any step"
                        )
                    if member in owner:
                        errors.append(
                            f"species.cycles.{label}: {member} is also in cycle {owner[member]}"
                        )
                    owner.setdefault(member, label)
                cycles[str(label)] = tuple(members)
    correction = species_block.get("vibrational_correction", "RRHO")
    if correction not in ("RRHO", "Grimme"):
        errors.append(
            f"species.vibrational_correction: must be RRHO or Grimme, got {correction!r}"
        )

    # Analyses
    ea = drc = mk = None
    raw = analyses.get("activation_energy")
    if raw is not None:
        where = "analyses.activation_energy"
        initial = _concentrations(raw.get("initial_M"), f"{where}.initial_M", errors)
        temps = _range(
            raw.get("temperatures_K"),
            f"{where}.temperatures_K",
            errors,
            minimum_points=3,
        )
        total = _number(
            raw.get("simulation_time_s"), f"{where}.simulation_time_s", errors
        )
        sampling = _number(
            raw.get("sampling_time_h"), f"{where}.sampling_time_h", errors
        )
        if None not in (total, sampling, output_step):
            _on_grid(sampling, output_step, total, f"{where}.sampling_time_h", errors)
        ea = ActivationEnergy(initial, temps or (), total, sampling)
    raw = analyses.get("degree_of_rate_control")
    if raw is not None:
        where = "analyses.degree_of_rate_control"
        if "initial_M" in raw:
            initial = _concentrations(raw["initial_M"], f"{where}.initial_M", errors)
        elif ea is not None:
            initial = dict(ea.initial_M)
        else:
            initial = {}
            errors.append(
                f"{where}.initial_M: required when activation_energy is not requested"
            )
        cores = _number(
            raw.get("cores", os.cpu_count() or 1),
            f"{where}.cores",
            errors,
            integer=True,
        )
        drc = DegreeOfRateControl(
            initial,
            _number(raw.get("e_shift_kcal_mol"), f"{where}.e_shift_kcal_mol", errors),
            _number(raw.get("simulation_time_s"), f"{where}.simulation_time_s", errors),
            cores,
        )
    raw = analyses.get("microkinetics")
    if raw is not None:
        where = "analyses.microkinetics"
        initial = _concentrations(raw.get("initial_M"), f"{where}.initial_M", errors)
        total = _number(
            raw.get("simulation_time_s"), f"{where}.simulation_time_s", errors
        )
        snapshot = _number(
            raw.get("catalyst_snapshot_h"), f"{where}.catalyst_snapshot_h", errors
        )
        conversion = _number(raw.get("conversion"), f"{where}.conversion", errors)
        if conversion is not None and not conversion < 1:
            errors.append(
                f"{where}.conversion: must be between 0 and 1, got {conversion}"
            )
        if None not in (total, snapshot, output_step):
            _on_grid(
                snapshot, output_step, total, f"{where}.catalyst_snapshot_h", errors
            )
        plot = raw.get("plot_species", [])
        if not isinstance(plot, list):
            errors.append(f"{where}.plot_species: must be a list of species")
            plot = []
        for name in plot:
            if names and name not in names:
                errors.append(
                    f"{where}.plot_species: {name} does not appear in any step"
                )
        if not cycles:
            errors.append("species.cycles: required for microkinetics")
        mk = Microkinetics(initial, total, snapshot, conversion, tuple(plot))

    # Names used in initial concentrations
    for label, analysis in (
        ("activation_energy", ea),
        ("degree_of_rate_control", drc),
        ("microkinetics", mk),
    ):
        if analysis is None or (
            label == "degree_of_rate_control" and "initial_M" not in analyses[label]
        ):
            continue
        where = f"analyses.{label}.initial_M"
        for name in analysis.initial_M:
            if names and name not in names:
                errors.append(f"{where}.{name}: does not appear in any step")
            if name == studied:
                errors.append(
                    f"{where}.{name}: the studied species is set by conditions.studied_range_M"
                )
            if name == product:
                errors.append(f"{where}.{name}: the product always starts at 0")
        members = {m for ms in cycles.values() for m in ms}
        if cycles and analysis.initial_M and not members & set(analysis.initial_M):
            errors.append(
                f"{where}: give the initial concentration of at least one catalyst intermediate"
            )

    # Energy sources
    needed = (
        sorted({temperature, *(ea.temperatures_K if ea else ())} - {None})
        if temperature
        else []
    )
    files = energies = None
    everything = names + [ts for ts in transition_states if ts not in names]
    if "files" in species_block and "energies_kcal_mol" not in species_block:
        files = (path.parent / str(species_block["files"])).resolve()
        if check_energy_sources:
            if not files.is_dir():
                errors.append(f"species.files: folder {files} not found")
            else:
                for name in everything:
                    if not (files / f"{name}.out").is_file():
                        errors.append(
                            f"no energy for {name}: {species_block['files']}/{name}.out not found"
                        )
            thermochange = os.environ.get("thermochange")
            if not thermochange:
                errors.append(
                    "$thermochange is not set: export thermochange=/path/to/thermochange"
                )
            elif not os.path.isfile(os.path.join(thermochange, THERMOCHANGE_SCRIPT)):
                errors.append(
                    f"$thermochange={thermochange} does not contain {THERMOCHANGE_SCRIPT}"
                )
    elif "energies_kcal_mol" in species_block and "files" not in species_block:
        energies = _energies(
            species_block["energies_kcal_mol"],
            "species.energies_kcal_mol",
            errors,
            needed,
            temperature,
            ea is not None,
        )
        for name in everything:
            if name not in energies:
                errors.append(
                    f"no energy for {name}: add it to species.energies_kcal_mol"
                )

    # Warnings
    if cycles:
        members = {m for ms in cycles.values() for m in ms}
        overall_names = set(overall.species()) if overall else set()
        loose = [
            n
            for n in names
            if n not in members and n not in overall_names and n != studied
        ]
        if loose:
            warnings.append(
                f"warning: {', '.join(loose)} belong to no cycle, so the catalyst distribution does not count them"
            )

    if errors:
        raise StudyError(errors)
    return Study(
        path=path,
        name=str(doc.get("name", path.stem)),
        files=files,
        energies_kcal_mol=energies,
        product=product,
        overall=overall,
        cycles=cycles,
        vibrational_correction=correction,
        steps=tuple(steps),
        temperature_K=temperature,
        studied_species=studied,
        studied_concentrations_M=concentrations,
        output_step_s=output_step,
        activation_energy=ea,
        degree_of_rate_control=drc,
        microkinetics=mk,
        warnings=tuple(warnings),
    )
```

- [ ] **Step 4: Run the tests**

Run: `python tests/test_study.py && python tests/test_steps.py && python tests/test_study_yaml.py`
Expected: `ok` three times. Each message asserted in `tests/test_study.py` is promised by spec §5: change the code, not the message, if one fails.

- [ ] **Step 5: Commit**

```bash
ruff check . && ruff format --check .
git add study.py tests/test_study.py
git commit -m "feat: load and check study files, reporting every problem with its location"
```

---

### Task 5: Gibbs energies and barrier tables

**Files:**
- Create: `thermochemistry.py`
- Create: `tests/test_thermochemistry.py`

**Interfaces:**
- Consumes: `Study` (Task 4); `GibbsEnergyCalculator`, `G_COMPOUNDS_OUTPUT_DIR_NAME`, `REACTION_DF_OUTPUT_DIR_NAME` (existing `calculating_G_for_microkinetics.py`); `AuxiliaryFunctions.compute_pressure_value` (existing).
- Produces:
  - `write_barrier_tables(study) -> list[str]`: writes `G_values_of_compounds/G_values_at_{T}K_{P:.5e}atm.csv` and `G_values_of_reactions/reaction_df_{T}K_{P:.5e}atm.csv` into the **current folder** for every `study.temperatures()`; keeps tables already present; returns warnings and the overall ΔG line.
  - `gibbs_energies(study, temperature_K) -> dict[str, float]`, `barrier_table(study, gibbs) -> pandas.DataFrame` (columns `Rx, TS, Gdir, Ginv`), `barrier_table_name(temperature_K) -> str`, `gibbs_from_file(out_file, temperature_K, correction, thermochange) -> float`, `class ThermochemistryError(RuntimeError)`.

- [ ] **Step 1: Write the failing tests** (`tests/test_thermochemistry.py`)

The fake thermochange writes `temp_summary.temp` into its working folder, as the real one does, and records its arguments, so the tests check the clean-up, the exact command line (including a folder name with a space) and the `-g` flag.

```python
"""Gibbs energies and barrier tables from typed energies and from (fake) thermochange.

python tests/test_thermochemistry.py
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from study import load_study
from thermochemistry import (
    ThermochemistryError,
    barrier_table_name,
    write_barrier_tables,
)

TYPED = {
    "species": {
        "energies_kcal_mol": {"S": 0.0, "P": -10.0, "C1": 0.0, "C2": -2.0, "TS1": 15.0},
        "product": "P",
        "overall_reaction": "S <=> P",
        "cycles": {"cat": ["C1", "C2"]},
    },
    "steps": ["C1 + S <=> C2", "C2 <=> C1 + P  via TS1"],
    "conditions": {
        "temperature_K": 300,
        "studied_species": "S",
        "studied_range_M": {"from": 0.01, "to": 0.1, "points": 3},
    },
    "analyses": {
        "microkinetics": {
            "initial_M": {"C1": 1e-3},
            "simulation_time_s": 1000,
            "catalyst_snapshot_h": 0.01,
            "conversion": 0.99,
            "plot_species": ["C1", "C2"],
        }
    },
}

FAKE_THERMOCHANGE = """#!/bin/bash
echo "$@" >> "$FAKE_LOG"
[ "$1" = "-g" ] && shift
echo leftover > temp_summary.temp
if [ "$FAKE_MODE" = zero ]; then
  printf "%s\\t-1.000000\\t-1.000000\\t0.000000\\n" "$1"
else
  printf "%s\\t-1.000000\\t-1.000000\\t-1.001000\\n" "$1"
fi
"""


def study_from(doc, files=()):
    folder = Path(tempfile.mkdtemp(prefix="study "))
    for name in files:
        (folder / name).parent.mkdir(parents=True, exist_ok=True)
        (folder / name).touch()
    (folder / "study.yaml").write_text(yaml.safe_dump(doc))
    return load_study(folder / "study.yaml")


def in_new_folder():
    os.chdir(tempfile.mkdtemp())


def test_typed_barriers_match_hand_calculation():
    study = study_from(TYPED)
    in_new_folder()
    messages = write_barrier_tables(study)
    table = pd.read_csv(Path("G_values_of_reactions") / barrier_table_name(300.0))
    # C1 + S <=> C2, barrierless: TS at max(0 + 0, -2) + 4 = 4
    # C2 <=> C1 + P via TS1: Gdir = 15 - (-2) = 17, Ginv = 15 - (0 - 10) = 25
    assert table["Rx"].tolist() == ["C1 + S = C2", "C2 = C1 + P"]
    assert table["TS"].tolist() == ["-", "TS1"]
    assert table["Gdir"].tolist() == [4.0, 17.0] and table["Ginv"].tolist() == [
        6.0,
        25.0,
    ], table
    assert "Overall reaction S <=> P: ΔG = -10.0 kcal/mol at 300.0 K" in messages, (
        messages
    )


def test_barrier_table_names_use_float_temperatures():
    study = study_from(TYPED)
    in_new_folder()
    write_barrier_tables(study)
    assert os.listdir("G_values_of_reactions") == [
        "reaction_df_300.0K_2.46172e+01atm.csv"
    ]


def test_coefficient_multiplies_the_energy():
    doc = {**TYPED, "steps": ["2 S <=> C2 via TS1", "C2 <=> C1 + P"]}
    doc["species"] = {
        **TYPED["species"],
        "overall_reaction": "2 S <=> P",
        "energies_kcal_mol": {"S": 1.0, "P": -10.0, "C1": 0.0, "C2": -3.0, "TS1": 10.0},
    }
    study = study_from(doc)
    in_new_folder()
    write_barrier_tables(study)
    table = pd.read_csv(Path("G_values_of_reactions") / barrier_table_name(300.0))
    assert table.loc[0, "Rx"] == "S + S = C2"
    assert table.loc[0, "Gdir"] == 10.0 - 2 * 1.0 and table.loc[0, "Ginv"] == 10.0 - (
        -3.0
    )


def test_negative_barrier_is_a_warning():
    doc = {
        **TYPED,
        "species": {
            **TYPED["species"],
            "energies_kcal_mol": {**TYPED["species"]["energies_kcal_mol"], "TS1": -5.0},
        },
    }
    study = study_from(doc)
    in_new_folder()
    messages = write_barrier_tables(study)
    assert any(
        "steps[1] (C2 <=> C1 + P via TS1) has a negative forward barrier (-3.0 kcal/mol)"
        in m
        for m in messages
    ), messages


def fake_thermochange(mode):
    root = Path(tempfile.mkdtemp())
    (root / "formatters").mkdir()
    script = root / "formatters" / "formatted_energy_outputter.sh"
    script.write_text(FAKE_THERMOCHANGE)
    os.environ.update(
        thermochange=str(root), FAKE_MODE=mode, FAKE_LOG=str(root / "calls.log")
    )
    return root


def files_study(correction="RRHO"):
    doc = {
        **TYPED,
        "species": {
            k: v for k, v in TYPED["species"].items() if k != "energies_kcal_mol"
        },
    }
    doc["species"].update(files="out", vibrational_correction=correction)
    return study_from(
        doc, files=[f"out/{n}.out" for n in ("S", "P", "C1", "C2", "TS1")]
    )


def test_file_energies_are_converted_and_leave_no_temporary_files():
    root = fake_thermochange("ok")
    study = files_study()
    in_new_folder()
    write_barrier_tables(study)
    compounds = pd.read_csv(
        Path("G_values_of_compounds") / "G_values_at_300.0K_2.46172e+01atm.csv",
        index_col=0,
    )
    assert abs(compounds.loc["S", "Gibbs Free Energies"] - (-1.001 * 627.509)) < 1e-9
    assert sorted(os.listdir(".")) == [
        "G_values_of_compounds",
        "G_values_of_reactions",
    ], os.listdir(".")
    calls = (root / "calls.log").read_text().splitlines()
    expected = f" 300.0 {0.082057366080960 * 300.0}"  # the same strings the old shell pipeline passed
    assert len(calls) == 5 and all(c.endswith(expected) for c in calls), calls
    assert all(
        " " in c.split(".out")[0] for c in calls
    )  # a folder name with a space stays one argument


def test_grimme_passes_the_g_flag():
    root = fake_thermochange("ok")
    study = files_study("Grimme")
    in_new_folder()
    write_barrier_tables(study)
    assert all(
        c.startswith("-g ") for c in (root / "calls.log").read_text().splitlines()
    )


def test_zero_energy_from_thermochange_stops_the_run():
    fake_thermochange("zero")
    study = files_study()
    in_new_folder()
    try:
        write_barrier_tables(study)
    except ThermochemistryError as error:
        assert "corrected G of 0.0 hartree" in str(error), str(error)
    else:
        raise AssertionError("expected ThermochemistryError")
    assert not Path("G_values_of_reactions").exists() or not os.listdir(
        "G_values_of_reactions"
    )


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
```

- [ ] **Step 2: Run to see it fail**

Run: `python tests/test_thermochemistry.py`
Expected: `ModuleNotFoundError: No module named 'thermochemistry'`

- [ ] **Step 3: Implement `thermochemistry.py`**

```python
"""Gibbs energies of a study's species and the barrier tables the COPASI simulations read"""

import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from auxiliary_functions import AuxiliaryFunctions
from calculating_G_for_microkinetics import (
    G_COMPOUNDS_OUTPUT_DIR_NAME,
    REACTION_DF_OUTPUT_DIR_NAME,
    GibbsEnergyCalculator,
)

HARTREE_TO_KCAL_MOL = 627.509
# A temperature/pressure correction moves G by millihartrees; more than this means thermochange failed
MAX_CORRECTION_HARTREE = 0.1


class ThermochemistryError(RuntimeError):
    """thermochange failed or gave an implausible Gibbs energy"""


def barrier_table_name(temperature_K):
    """File name the simulations look up for the barriers at temperature_K"""
    pressure = AuxiliaryFunctions.compute_pressure_value(temperature_K)
    return f"reaction_df_{temperature_K}K_{pressure:.5e}atm.csv"


def compounds_table_name(temperature_K):
    pressure = AuxiliaryFunctions.compute_pressure_value(temperature_K)
    return f"G_values_at_{temperature_K}K_{pressure:.5e}atm.csv"


def gibbs_from_file(out_file, temperature_K, correction, thermochange):
    """G (kcal/mol, 1 M) of one output file at temperature_K, through thermochange"""
    pressure = AuxiliaryFunctions.compute_pressure_value(temperature_K)
    script = Path(thermochange) / "formatters" / "formatted_energy_outputter.sh"
    flags = ["-g"] if correction == "Grimme" else []
    command = [
        "bash",
        str(script),
        *flags,
        str(out_file),
        str(temperature_K),
        str(pressure),
    ]
    # thermochange writes temp_summary.temp to the current folder and does not always remove it
    with tempfile.TemporaryDirectory() as scratch:
        result = subprocess.run(
            command, cwd=scratch, capture_output=True, text=True, check=False
        )
    name = Path(out_file).stem
    lines = result.stdout.strip().splitlines()
    fields = lines[-1].split("\t") if lines else []
    try:
        file_g, corrected = float(fields[2]), float(fields[3])
    except (IndexError, ValueError) as error:
        raise ThermochemistryError(
            f"thermochange gave no Gibbs energy for {name}:\n{result.stdout}{result.stderr}"
        ) from error
    if corrected == 0.0 or abs(corrected - file_g) > MAX_CORRECTION_HARTREE:
        raise ThermochemistryError(
            f"thermochange gave a corrected G of {corrected} hartree for {name} "
            f"(the file's own G is {file_g}); its error output:\n{result.stderr}"
        )
    return corrected * HARTREE_TO_KCAL_MOL


def _typed(study, name, temperature_K):
    for known, value in study.energies_kcal_mol[name].items():
        if abs(known - temperature_K) < 1e-6:
            return value
    raise ThermochemistryError(f"no typed energy for {name} at {temperature_K} K")


def gibbs_energies(study, temperature_K):
    """{species or transition state: G in kcal/mol} at temperature_K"""
    names = study.species() + [
        ts for ts in study.transition_states() if ts not in study.species()
    ]
    if study.files is None:
        return {name: _typed(study, name, temperature_K) for name in names}
    thermochange = os.environ["thermochange"]
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 1) as pool:
        values = list(
            pool.map(
                lambda name: gibbs_from_file(
                    study.files / f"{name}.out",
                    temperature_K,
                    study.vibrational_correction,
                    thermochange,
                ),
                names,
            )
        )
    return dict(zip(names, values))


def barrier_table(study, gibbs):
    """Rx, TS, Gdir, Ginv for every step, in the format copasi_helper reads"""
    table = pd.DataFrame(
        {"Rx": study.reactions(), "TS": [step.ts or "-" for step in study.steps]}
    )
    energies = pd.DataFrame({"Gibbs Free Energies": pd.Series(gibbs)})
    gdir, ginv = GibbsEnergyCalculator(
        energies
    ).calculate_direct_inverse_reactions_gibbs_free_energies(table)
    table["Gdir"], table["Ginv"] = gdir, ginv
    return table


def negative_barrier_warnings(study, table):
    messages = []
    for i, (step, gdir, ginv) in enumerate(
        zip(study.steps, table["Gdir"], table["Ginv"])
    ):
        for value, direction in ((gdir, "forward"), (ginv, "reverse")):
            if value < 0:
                messages.append(
                    f"warning: steps[{i}] ({step.text}) has a negative {direction} barrier ({value:.1f} kcal/mol)"
                )
    return messages


def overall_reaction_energy(study, gibbs):
    """ΔG of the overall reaction from the species energies (kcal/mol)"""
    products = sum(n * gibbs[name] for name, n in study.overall.products)
    reactants = sum(n * gibbs[name] for name, n in study.overall.reactants)
    return products - reactants


def write_barrier_tables(study):
    """Writes the compound and barrier tables for every needed temperature into the current folder.

    Tables already present are kept (the results folder is per study, guarded by run_info.json).
    Returns the messages to show: warnings and the overall reaction energy at the working temperature.
    """
    messages = []
    for temperature in study.temperatures():
        reactions_path = Path(REACTION_DF_OUTPUT_DIR_NAME) / barrier_table_name(
            temperature
        )
        compounds_path = Path(G_COMPOUNDS_OUTPUT_DIR_NAME) / compounds_table_name(
            temperature
        )
        if reactions_path.is_file() and compounds_path.is_file():
            continue
        print(f"Computing Gibbs energies at {temperature} K...")
        gibbs = gibbs_energies(study, temperature)
        table = barrier_table(study, gibbs)
        compounds_path.parent.mkdir(exist_ok=True)
        reactions_path.parent.mkdir(exist_ok=True)
        pd.DataFrame({"Gibbs Free Energies": pd.Series(gibbs)}).rename_axis(
            "Compounds"
        ).to_csv(compounds_path)
        partial = reactions_path.with_suffix(".partial")
        table.to_csv(partial, index=False)
        partial.replace(
            reactions_path
        )  # written last and atomically: its presence means complete
        messages += negative_barrier_warnings(study, table)
        if temperature == study.temperature_K:
            delta = overall_reaction_energy(study, gibbs)
            messages.append(
                f"Overall reaction {study.overall.text}: ΔG = {delta:.1f} kcal/mol at {temperature} K"
            )
    return messages
```

- [ ] **Step 4: Run the tests**

Run: `python tests/test_thermochemistry.py`
Expected: `ok`

- [ ] **Step 5: Commit**

```bash
ruff check . && ruff format --check .
git add thermochemistry.py tests/test_thermochemistry.py
git commit -m "feat: Gibbs energies from thermochange or typed values, and barrier tables per temperature"
```

---

### Task 6: Generalize the analyses

Removes the assumptions that block a second system: the G step inside the analyses, `CURRENT_DIRECTORY` fixed at import, the product named `prod`, the cycle regex, the 1:1 yield, and fixed figure grids. `main.py` keeps working (it now writes the barrier tables up front) until Task 8 replaces it.

**Files:**
- Modify: `file_operations.py`, `apparent_activation_energy.py`, `microkinetics_simulation.py`, `auxiliary_functions.py`, `calculating_G_for_microkinetics.py`, `plotting_functions.py`, `main.py` (patch below)
- Create: `tests/test_generalized_pipeline.py`

**Interfaces:**
- Produces:
  - `file_operations.run_directory() -> str` (the current folder; replaces `CURRENT_DIRECTORY` everywhere)
  - `DRCAnalysis(..., cores=1, product="prod", time_step=1)`
  - `MicroKinetics(temperature_value, total_simulation_time, reactant_concentration_array, reactant_to_study, c0, cycles: dict[str, list[str]], max_product_M: list[float], compounds_to_plot, catalyst_concentration_time, percentage_of_convertion, time_step, product="prod")`
  - `AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(simulation_dfs, reactant_concentration_array, reactant_to_study, thresholds: list[float], product="prod")`
  - `PlotFunctions.grid(panels, max_columns=5) -> (nrows, ncols)`
  - `ApparentEaAnalysis` no longer computes Gibbs energies: the barrier tables must exist before it runs.

- [ ] **Step 1: Write the failing tests** (`tests/test_generalized_pipeline.py`)

```python
"""The analyses no longer assume a product named prod, regex-named cycles or a 1:1 yield.

Uses a fake copasi_parser, so it runs without COPASI:
    python tests/test_generalized_pipeline.py
"""

import os
import sys
import tempfile
import types

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.modules.setdefault("copasi_parser", types.ModuleType("copasi_parser"))

from auxiliary_functions import AuxiliaryFunctions
from plotting_functions import PlotFunctions


def test_conversion_threshold_is_per_concentration():
    time = np.arange(0, 11) / 10
    sims = [
        pd.DataFrame({"time": time, "P": time * 1.0}),
        pd.DataFrame({"time": time, "P": time * 2.0}),
    ]
    # First simulation must pass 0.5, second must pass 0.5 too: thresholds differ per concentration
    times = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
        sims, [0.01, 0.1], "S", [0.45, 1.5], "P"
    )
    assert times == [(0.5, 0.01), (0.8, 0.1)], times


def test_grid_sizes():
    assert PlotFunctions.grid(1) == (1, 1)
    assert PlotFunctions.grid(5) == (1, 5)
    assert PlotFunctions.grid(6) == (2, 5)
    assert PlotFunctions.grid(22) == (5, 5)
    assert PlotFunctions.grid(0) == (1, 1)


def test_results_follow_the_current_folder():
    from file_operations import FileOperations, run_directory

    folder = tempfile.mkdtemp()
    os.chdir(folder)
    assert run_directory() == os.getcwd()
    open("table.csv", "w").close()
    FileOperations.move_to_output_directory("out", "table.csv")
    assert os.path.isfile(os.path.join(folder, "out", "table.csv"))


def test_drc_uses_the_product_and_time_step():
    calls = []

    class FakeModel:
        def self_destruct(self):
            pass

    cpx = sys.modules["copasi_parser"]
    cpx.prepare_copasi_model = lambda **kwargs: FakeModel()
    cpx.drc_calc = lambda base_model, **kwargs: (
        calls.append(kwargs) or (np.array([1.0]), 1.0)
    )
    os.chdir(tempfile.mkdtemp())
    from apparent_activation_energy import DRCAnalysis

    DRCAnalysis(
        np.array([1e-3]),
        np.array([300.0]),
        100,
        {"C1": 1e-3},
        "S",
        ["A = B"],
        0.1,
        1,
        product="P",
        time_step=5,
    ).calculate_degree_of_rate_control()
    assert {c["target_spc"] for c in calls} == {"P"} and {
        c["time_step"] for c in calls
    } == {5}, calls


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
```

- [ ] **Step 2: Run to see it fail**

Run: `python tests/test_generalized_pipeline.py`
Expected: `TypeError` (unexpected argument `product`, or `thresholds` not subscriptable) or `AttributeError: type object 'PlotFunctions' has no attribute 'grid'`

- [ ] **Step 3: Apply the patch**

Save the block below exactly as `task6.patch` in the repository root, then apply it (`git apply --check` first):

```diff
--- a/file_operations.py
+++ b/file_operations.py
@@ -2,7 +2,10 @@
 
 import os
 
-CURRENT_DIRECTORY = os.getcwd()
+
+def run_directory():
+    """Folder all results are written to: the current folder (microkatc.py runs inside results/)"""
+    return os.getcwd()
 
 
 class FileOperations:
@@ -11,8 +14,8 @@
     @staticmethod
     def move_to_output_directory(output_dir_name, file_name):
         """Moves a file to the specified directory"""
-        current_file_path = os.path.join(CURRENT_DIRECTORY, file_name)
-        desired_file_path = os.path.join(CURRENT_DIRECTORY, output_dir_name, file_name)
+        current_file_path = os.path.join(run_directory(), file_name)
+        desired_file_path = os.path.join(run_directory(), output_dir_name, file_name)
 
         os.makedirs(os.path.dirname(desired_file_path), exist_ok=True)
 
--- a/apparent_activation_energy.py
+++ b/apparent_activation_energy.py
@@ -10,7 +10,7 @@
 
 from auxiliary_functions import AuxiliaryFunctions
 from calculating_G_for_microkinetics import REACTION_DF_OUTPUT_DIR_NAME
-from file_operations import CURRENT_DIRECTORY, FileOperations
+from file_operations import FileOperations, run_directory
 from microkinetics_simulation import SIMULATIONS_OUTPUT_DIR_NAME, SimulationHandler
 from plotting_functions import PlotFunctions
 
@@ -40,8 +40,8 @@
             c0_copy[reactant_to_study] = initial_concentration
 
             for T_value in T_values_array:
+                # The barrier tables are written beforehand (thermochemistry.write_barrier_tables)
                 pressure_value = AuxiliaryFunctions.compute_pressure_value(T_value)
-                AuxiliaryFunctions.calculate_G_values(T_value, pressure_value)
 
                 simulation_calculator = SimulationHandler(
                     T_value,
@@ -256,7 +256,7 @@
 
         try:
             df_flux_path = os.path.join(
-                CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, self.df_flux_filename
+                run_directory(), SIMULATIONS_OUTPUT_DIR_NAME, self.df_flux_filename
             )
             self._df_flux = pd.read_csv(df_flux_path)
             print(
@@ -268,7 +268,7 @@
 
         try:
             df_rate_path = os.path.join(
-                CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, self.df_rate_filename
+                run_directory(), SIMULATIONS_OUTPUT_DIR_NAME, self.df_rate_filename
             )
             self._df_rate = pd.read_csv(df_rate_path)
             print(
@@ -401,6 +401,8 @@
         reactions,
         e_shift,
         cores=1,
+        product="prod",
+        time_step=1,
     ):
         self.reactant_concentration_array = reactant_concentration_array
         self.T_values_array = T_values_array
@@ -410,6 +412,8 @@
         self.reactions = reactions
         self.e_shift = e_shift
         self.cores = cores
+        self.product = product
+        self.time_step = time_step
 
         self.plot_function = PlotFunctions()
 
@@ -422,12 +426,14 @@
             self.reactions,
             self.e_shift,
             "central difference",
+            self.product,
+            self.time_step,
         )
         self.df_drc_filename = f"df_drc_T_range_{self.T_values_array[0]}K_{self.T_values_array[-1]}K_C({self.reactant_to_study})_{self.reactant_concentration_array[0]}M_{self.reactant_concentration_array[-1]}M_{inputs_hash}.csv"
 
         try:
             df_drc_path = os.path.join(
-                CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, self.df_drc_filename
+                run_directory(), SIMULATIONS_OUTPUT_DIR_NAME, self.df_drc_filename
             )
             self._df_drc = pd.read_csv(df_drc_path)
             print(
@@ -457,7 +463,7 @@
         for T_value in self.T_values_array:
             pressure_value = AuxiliaryFunctions.compute_pressure_value(T_value)
             datafile = os.path.join(
-                CURRENT_DIRECTORY,
+                run_directory(),
                 REACTION_DF_OUTPUT_DIR_NAME,
                 f"reaction_df_{T_value}K_{pressure_value:.5e}atm.csv",
             )
@@ -515,9 +521,9 @@
             temp=T_value,
             c0=c0,
             total_time=self.total_simulation_time,
-            time_step=1,
+            time_step=self.time_step,
             target_time=self.total_simulation_time / 2,
-            target_spc="prod",
+            target_spc=self.product,
             e_shift=shift,
             cores=self.cores,
         )
--- a/microkinetics_simulation.py
+++ b/microkinetics_simulation.py
@@ -7,7 +7,7 @@
 
 from auxiliary_functions import AuxiliaryFunctions
 from calculating_G_for_microkinetics import REACTION_DF_OUTPUT_DIR_NAME
-from file_operations import CURRENT_DIRECTORY, FileOperations
+from file_operations import FileOperations, run_directory
 from plotting_functions import PlotFunctions
 
 CONVERT_SECONDS_TO_HOURS = 1 / 3600
@@ -24,7 +24,7 @@
         self.total_simulation_time = total_simulation_time
         self.time_step = time_step
         self.datafile = os.path.join(
-            CURRENT_DIRECTORY,
+            run_directory(),
             REACTION_DF_OUTPUT_DIR_NAME,
             f"reaction_df_{self.temperature_value}K_{self.pressure_value:.5e}atm.csv",
         )
@@ -91,7 +91,7 @@
         )
 
         simulation_file_path = os.path.join(
-            CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, simulation_filename
+            run_directory(), SIMULATIONS_OUTPUT_DIR_NAME, simulation_filename
         )
 
         if os.path.exists(simulation_file_path):
@@ -124,12 +124,15 @@
         reactant_to_study,
         c0,
         cycles,
-        main_reactants,
+        max_product_M,
         compounds_to_plot,
         catalyst_concentration_time,
         percentage_of_convertion,
         time_step,
+        product="prod",
     ):
+        """cycles: {label: [intermediates]}; max_product_M: the largest product concentration the
+        overall reaction allows, one value per studied concentration (study.Study.max_product_M)"""
 
         self.temperature_value = temperature_value
         self.pressure_value = AuxiliaryFunctions.compute_pressure_value(
@@ -139,7 +142,8 @@
         self.total_simulation_time = total_simulation_time
         self.c0 = c0
         self.cycles = cycles
-        self.main_reactants = main_reactants
+        self.max_product_M = list(max_product_M)
+        self.product = product
         self.compounds_to_plot = compounds_to_plot
         self.reactant_concentration_array = reactant_concentration_array
         self.catalyst_concentration_time = catalyst_concentration_time
@@ -197,21 +201,18 @@
     def plot_catalyst_concentration_vs_reactant(self, log_x, figsize):
         """Plots and saves catalyst concentration in each cycle vs. initial concentration of the studied reactant"""
         fig_name = f"C(catalyst)_C({self.reactant_to_study})_{self.temperature_value}K_{self.pressure_value:.5e}atm.svg"
-        cycles_intermediates_dict = AuxiliaryFunctions.find_intermediates_of_cycle(
-            self.cycles
-        )
         catalyst_concentrations_per_cycle = (
             AuxiliaryFunctions.get_concentrations_of_catalyst(
                 self.simulations_dfs,
                 self.catalyst_concentration_time,
-                cycles_intermediates_dict,
+                self.cycles,
             )
         )
 
         self.plot_functions.plot_concentration_of_catalyst_versus_studied_reactant(
             self.reactant_concentration_array,
             catalyst_concentrations_per_cycle,
-            self.cycles,
+            list(self.cycles),
             self.reactant_to_study,
             figsize,
             fig_name,
@@ -222,18 +223,13 @@
     def plot_reactant_vs_product_conversion(self, log_x, figsize):
         """Plots and saves the time to reach the product conversion threshold vs. initial concentration of the studied reactant"""
         fig_name = f"C({self.reactant_to_study})_total_t_product_conversion_{self.temperature_value}K_{self.pressure_value:.5e}atm.svg"
-        product_conversion_threshold_concentration = (
-            AuxiliaryFunctions.find_limiting_reactant_concentration(
-                self.c0, self.main_reactants
-            )
-            * self.percentage_of_convertion
-        )
-
+        thresholds = [m * self.percentage_of_convertion for m in self.max_product_M]
         times_conv_reac_conc = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
             self.simulations_dfs,
             self.reactant_concentration_array,
             self.reactant_to_study,
-            product_conversion_threshold_concentration,
+            thresholds,
+            self.product,
         )
 
         self.plot_functions.plot_reactant_concentration_vs_product_conversion(
--- a/auxiliary_functions.py
+++ b/auxiliary_functions.py
@@ -12,7 +12,7 @@
     G_COMPOUNDS_OUTPUT_DIR_NAME,
     REACTION_DF_OUTPUT_DIR_NAME,
 )
-from file_operations import CURRENT_DIRECTORY
+from file_operations import run_directory
 
 R_L_atm_per_mol_K = 0.082057366080960
 
@@ -92,17 +92,19 @@
         simulation_dfs,
         reactant_concentration_array,
         reactant_to_study,
-        product_conversion_threshold_concentration,
+        thresholds,
+        product="prod",
     ):
-        """Returns the first times when concentration of a product is greater than a product concentration threshold"""
+        """First time each simulation's product concentration exceeds its threshold (one per simulation)"""
         times_conv_reac_conc = []
 
         for i, simulation_df in enumerate(simulation_dfs):
             try:
                 # Find the first time value when concentration of a product is greater than a set product concentration threshold
-                filtered_df = simulation_df.query(
-                    "prod > @product_conversion_threshold_concentration"
-                )
+                product_conversion_threshold_concentration = thresholds[i]
+                filtered_df = simulation_df[
+                    simulation_df[product] > product_conversion_threshold_concentration
+                ]
 
                 # First time at when product reach the product conversion threshold concentration
                 first_convergence_time = filtered_df.iloc[0]["time"]
@@ -133,17 +135,17 @@
         G_compounds_file_name = f"G_values_at_{temperature}K_{pressure:.5e}atm.csv"
 
         G_compounds_file_path = os.path.join(
-            CURRENT_DIRECTORY, G_COMPOUNDS_OUTPUT_DIR_NAME, G_compounds_file_name
+            run_directory(), G_COMPOUNDS_OUTPUT_DIR_NAME, G_compounds_file_name
         )
         reaction_df_file_path = os.path.join(
-            CURRENT_DIRECTORY, REACTION_DF_OUTPUT_DIR_NAME, reaction_df_file_name
+            run_directory(), REACTION_DF_OUTPUT_DIR_NAME, reaction_df_file_name
         )
         outputs = [G_compounds_file_path, reaction_df_file_path]
 
         # Calculations will be executed if either file for the specified T and P is missing
         if not all(os.path.exists(path) for path in outputs):
             # Path to the Bash script
-            script_path = os.path.join(CURRENT_DIRECTORY, "get_G_compounds.sh")
+            script_path = os.path.join(run_directory(), "get_G_compounds.sh")
             print(f"Creating and saving: {reaction_df_file_name}")
             print(f"Creating and saving: {G_compounds_file_name}")
             result = subprocess.run(
--- a/calculating_G_for_microkinetics.py
+++ b/calculating_G_for_microkinetics.py
@@ -5,7 +5,7 @@
 import pandas as pd
 
 from bash_parsing import BashParser
-from file_operations import CURRENT_DIRECTORY, FileOperations
+from file_operations import FileOperations, run_directory
 
 DIFFUSION_BARRIER = 4
 G_COMPOUNDS_OUTPUT_DIR_NAME = "G_values_of_compounds"
@@ -108,7 +108,7 @@
 
     def create_compounds_df(self, compounds):
         """Creates and saves dataframe with compounds in the system if it doesn't exist"""
-        compounds_df_path = os.path.join(CURRENT_DIRECTORY, self.df_compounds_file_name)
+        compounds_df_path = os.path.join(run_directory(), self.df_compounds_file_name)
         if not os.path.exists(compounds_df_path):
             print("Creating compounds.csv file...")
             compounds.to_csv(f"{self.df_compounds_file_name}", index=False)
--- a/plotting_functions.py
+++ b/plotting_functions.py
@@ -17,6 +17,12 @@
 
     def __init__(self):
         os.makedirs(IMAGES_SIMULATIONS_OUTPUT_DIR_NAME, exist_ok=True)
+
+    @staticmethod
+    def grid(panels, max_columns=5):
+        """(nrows, ncols) for a figure with the given number of panels"""
+        ncols = max(1, min(max_columns, panels))
+        return -(-max(panels, 1) // ncols), ncols
 
     @staticmethod
     def hide_unused_subplots(fig, axes):
--- a/main.py
+++ b/main.py
@@ -55,6 +55,12 @@
 
     # Specify time step for the microkinetics simulations
     time_step = 1
+
+    # Barrier tables for every temperature, before any simulation (the analyses no longer compute them)
+    for T_value in sorted({*T_values_array_Ea, temperature_value}):
+        AuxiliaryFunctions.calculate_G_values(
+            T_value, AuxiliaryFunctions.compute_pressure_value(T_value)
+        )
 
     # Perform MicroKinetics Analysis
     analysis1 = ApparentEaAnalysis(
@@ -179,8 +185,8 @@
         reactant_concentration_array,
         reactant_to_study,
         c0_2,
-        cycles,
-        main_reactants,
+        AuxiliaryFunctions.find_intermediates_of_cycle(cycles),
+        [min(c0_2[r] for r in main_reactants)] * len(reactant_concentration_array),
         compounds_to_plot,
         catalyst_concentration_time,
         percentage_of_convertion,
```

```bash
git apply --check task6.patch && git apply task6.patch && rm task6.patch
grep -rn "CURRENT_DIRECTORY" --include=*.py . || echo "no CURRENT_DIRECTORY left"
```
Expected: `no CURRENT_DIRECTORY left`. If `git apply --check` fails, the base is not `main` with #25–#27 merged and formatted by ruff 0.16.10: stop and fix the base.

- [ ] **Step 4: Run every test**

```bash
for t in tests/test_*.py; do echo "$t: $(python $t 2>&1 | tail -1)"; done
```
Expected: `ok` for every file (`test_paper_barriers.py` prints `skipped: ...` without `$thermochange`).

- [ ] **Step 5: Commit**

```bash
ruff check . && ruff format --check .
git add -A
git commit -m "refactor: analyses take the product, cycles, maximum yield and output step as inputs"
```

---

### Task 7: The `microkatc.py` command and the typed-energy example

**Files:**
- Create: `microkatc.py`
- Create: `examples/typed_energies/study.yaml`
- Create: `tests/test_microkatc.py`

**Interfaces:**
- Consumes: `load_study`, `Study` (Task 4); `write_barrier_tables` (Task 5); the Task 6 signatures.
- Produces:
  - `main(argv: list[str] | None) -> int` — `run <study> [--fresh]` and `check <study>`; exit codes 0/1/2.
  - `prepare_results(study, fresh=False) -> (Path, dict)`, `run_analyses(study)` (runs in the current folder), `initial_concentrations(study, initial_M, studied_value) -> dict` (used by `readme_figures.py` and `tests/check_reproduction.py` so their simulations hit the same cache), `study_hash(study) -> str`, `class StaleResultsError(Exception)`.
  - `results/` contents: `run_info.json`, `microkinetics.json` (same format as the Task 1 baseline), `G_values_of_*`, `microkinetics_simulations*`, `model.cps`.

**Measured constraint:** COPASI resolves relative output paths against the first folder it saved a model in during a process, so a second study in the same Python process writes into the first study's `results/` (reproduced 2026-10-02). The tests therefore run each study through a subprocess, as users do; the docstring of `microkatc.py` states the rule.

- [ ] **Step 1: Create the typed-energy example** (`examples/typed_energies/study.yaml`)

```yaml
# MicroKatc study with typed energies: one catalytic cycle, no output files needed.
# C1 binds the substrate S (barrierless) and releases the product P over TS1.
name: Typed-energy example, one catalytic cycle

species:
  energies_kcal_mol:               # Gibbs energies at 300 K and 1 M, kcal/mol, any common reference
    S: 0.0
    P: -10.0
    C1: 0.0
    C2: -2.0
    TS1: 15.0
  product: P
  overall_reaction: S <=> P
  cycles:
    cat: [C1, C2]

steps:
  - C1 + S <=> C2                  # barrierless: 4 kcal/mol above the higher side
  - C2 <=> C1 + P  via TS1

conditions:
  temperature_K: 300
  studied_species: S
  studied_range_M: {from: 0.01, to: 0.1, points: 3}

analyses:
  microkinetics:
    initial_M: {C1: 1e-3}
    simulation_time_s: 1000
    catalyst_snapshot_h: 0.01      # 36 s
    conversion: 0.99
    plot_species: [C1, C2, S, P]
```

Barriers by hand (used in Task 5's test as well): `C1 + S <=> C2` barrierless, TS at max(0 + 0, −2) + 4 = 4, so Gdir = 4, Ginv = 6; `C2 <=> C1 + P via TS1`: Gdir = 15 − (−2) = 17, Ginv = 15 − (0 − 10) = 25.

- [ ] **Step 2: Write the failing tests** (`tests/test_microkatc.py`)

```python
"""The microkatc.py command on the typed-energy example (real COPASI, runs in a few seconds).

python tests/test_microkatc.py
"""

import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import microkatc

EXAMPLE = ROOT / "examples" / "typed_energies" / "study.yaml"


def copy_example():
    folder = Path(
        tempfile.mkdtemp(prefix="my study ")
    )  # a space in the path on purpose
    shutil.copy(EXAMPLE, folder / "study.yaml")
    return folder / "study.yaml"


def run(*argv):
    """microkatc.py in its own process, as users run it (one study per process, see microkatc.py)"""
    result = subprocess.run(
        [sys.executable, str(ROOT / "microkatc.py"), *argv],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.returncode, result.stdout, result.stderr


def run_here(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = microkatc.main(list(argv))
    return code, out.getvalue(), err.getvalue()


def test_check_reports_ok_and_errors():
    study = copy_example()
    code, out, _ = run("check", str(study))
    assert code == 0 and "OK (2 steps, 4 species)" in out, out
    study.write_text(study.read_text().replace("product: P", "product: Q"))
    code, _, err = run("check", str(study))
    assert code == 1 and "species.product: Q does not appear in any step" in err, err
    assert not (study.parent / "results").exists()  # check never writes results


def test_run_from_other_folder():
    study = copy_example()
    started_in = tempfile.mkdtemp()
    result = subprocess.run(
        [sys.executable, str(ROOT / "microkatc.py"), "run", str(study)],
        cwd=started_in,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    results = study.parent / "results"
    assert sorted(p.name for p in results.iterdir()) == [
        "G_values_of_compounds",
        "G_values_of_reactions",
        "microkinetics.json",
        "microkinetics_simulations",
        "microkinetics_simulations_images",
        "model.cps",  # COPASI model of the last simulation, as before
        "run_info.json",
    ], sorted(p.name for p in results.iterdir())
    assert os.listdir(started_in) == []  # nothing written where the command was started
    summary = json.loads((results / "microkinetics.json").read_text())
    assert len(summary["t99_h"]) == 3 and all(
        abs(sum(x) - 1e-3) < 1e-9 for x in zip(*summary["catalyst_M"].values())
    )
    info = json.loads((results / "run_info.json").read_text())
    for key in (
        "study_sha256",
        "microkatc_commit",
        "python",
        "python-copasi",
        "numpy",
        "scipy",
        "started",
        "finished",
        "duration_s",
    ):
        assert key in info, key
    assert info["finished"] is not None


def test_rerun_same_study_is_allowed():
    study = copy_example()
    assert run("run", str(study))[0] == 0
    assert run("run", str(study))[0] == 0


def test_edited_study_stops_unless_fresh():
    study = copy_example()
    assert run("run", str(study))[0] == 0
    study.write_text(study.read_text().replace("TS1: 15.0", "TS1: 16.0"))
    code, _, err = run("run", str(study))
    assert (
        code == 1
        and "was made from a different version of study.yaml; delete it or run with --fresh"
        in err
    ), err
    assert run("run", str(study), "--fresh")[0] == 0


def test_failure_during_the_run_exits_2():
    study = copy_example()
    real = microkatc.run_analyses

    def broken(_study):
        raise RuntimeError("COPASI exploded")

    microkatc.run_analyses = broken
    try:
        code, _, err = run_here("run", str(study))
    finally:
        microkatc.run_analyses = real
    assert code == 2 and "run failed: COPASI exploded" in err, err


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
```

- [ ] **Step 3: Run to see it fail**

Run: `python tests/test_microkatc.py`
Expected: `ModuleNotFoundError: No module named 'microkatc'`

- [ ] **Step 4: Implement `microkatc.py`**

```python
"""MicroKatc command line.

    python microkatc.py run <study.yaml> [--fresh]    run every analysis the study requests
    python microkatc.py check <study.yaml>            check the study without running anything

Results go to a results/ folder next to the study file. Exit codes: 0 success, 1 invalid study
or stale results, 2 failure during the run.

Run one study per process: COPASI resolves relative output paths against the first folder it
saved a model in, so a second study run in the same Python process would write into the first
study's results/ (measured 2026-10-02).
"""

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from study import load_study
from study_yaml import StudyError

REPO = Path(__file__).resolve().parent
RUN_INFO = "run_info.json"


class StaleResultsError(Exception):
    """results/ was made from a different version of the study"""


def study_hash(study):
    """SHA-256 of the study file and of every energy input it uses"""
    digest = hashlib.sha256(study.path.read_bytes())
    if study.files is not None:
        for name in sorted(set(study.species()) | set(study.transition_states())):
            digest.update(name.encode())
            digest.update((study.files / f"{name}.out").read_bytes())
    return digest.hexdigest()


def _git_commit(folder):
    result = subprocess.run(
        ["git", "-C", str(folder), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() or None


def _version(package):
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def _write_run_info(results, info):
    (results / RUN_INFO).write_text(json.dumps(info, indent=2) + "\n")


def prepare_results(study, fresh=False):
    """Creates results/, or checks that the existing one was made from this exact study"""
    results = study.results_dir
    if fresh and results.exists():
        shutil.rmtree(results)
    digest = study_hash(study)
    info_path = results / RUN_INFO
    if info_path.is_file():
        if json.loads(info_path.read_text()).get("study_sha256") != digest:
            raise StaleResultsError(
                f"{results} was made from a different version of {study.path.name}; delete it or run with --fresh"
            )
    elif results.exists() and any(results.iterdir()):
        raise StaleResultsError(
            f"{results} has no {RUN_INFO}; delete it or run with --fresh"
        )
    results.mkdir(parents=True, exist_ok=True)
    thermochange = os.environ.get("thermochange")
    info = {
        "study": str(study.path),
        "study_sha256": digest,
        "microkatc_commit": _git_commit(REPO),
        "thermochange_commit": _git_commit(thermochange)
        if thermochange and study.files
        else None,
        "python": platform.python_version(),
        "python-copasi": _version("python-copasi"),
        "numpy": _version("numpy"),
        "scipy": _version("scipy"),
        "started": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "finished": None,
        "duration_s": None,
    }
    _write_run_info(results, info)
    return results, info


def initial_concentrations(study, initial_M, studied_value):
    """Initial concentrations for one simulation: the analysis' values, product at 0, studied species set"""
    return {
        **initial_M,
        study.product: 0.0,
        study.studied_species: float(studied_value),
    }


def run_analyses(study):
    """Runs every requested analysis in the current folder (results/)"""
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    from apparent_activation_energy import ApparentEaAnalysis, DRCAnalysis
    from auxiliary_functions import AuxiliaryFunctions
    from microkinetics_simulation import MicroKinetics
    from plotting_functions import PlotFunctions
    from thermochemistry import write_barrier_tables

    for message in write_barrier_tables(study):
        print(message)
    reactions = study.reactions()
    concentrations = np.array(study.studied_concentrations_M)
    T, studied, grid = study.temperature_K, study.studied_species, PlotFunctions.grid

    ea = study.activation_energy
    if ea is not None:
        analysis = ApparentEaAnalysis(
            T,
            np.array(ea.temperatures_K),
            concentrations,
            studied,
            initial_concentrations(study, ea.initial_M, concentrations[0]),
            ea.simulation_time_s,
            ea.sampling_time_h,
            reactions,
            study.output_step_s,
        )
        df_flux, df_rate = analysis.df_flux, analysis.df_rate
        analysis.plot_ln_ri_vs_1_over_T(
            df_flux, concentrations[0], *grid(len(reactions)), (15, 15)
        )
        forming = [
            r for r, s in zip(reactions, study.steps) if study.product in s.species()
        ]
        analysis.plot_Ea_vs_c0_flux_based(
            df_flux, forming, *grid(len(forming)), (15, 10), log_x=True
        )
        analysis.plot_Ea_vs_c0_rate_based(
            df_rate, [study.product], 1, 1, (15, 10), log_x=True
        )

    drc = study.degree_of_rate_control
    if drc is not None:
        analysis = DRCAnalysis(
            concentrations,
            np.array([T]),
            drc.simulation_time_s,
            initial_concentrations(study, drc.initial_M, concentrations[0]),
            studied,
            reactions,
            drc.e_shift_kcal_mol,
            drc.cores,
            product=study.product,
            time_step=study.output_step_s,
        )
        analysis.plot_drc_vs_c0(
            *grid(len(reactions)), (15, 15), analysis.df_drc, T, log_x=True
        )

    mk = study.microkinetics
    if mk is not None:
        cycles = {label: list(members) for label, members in study.cycles.items()}
        analysis = MicroKinetics(
            T,
            mk.simulation_time_s,
            concentrations,
            studied,
            initial_concentrations(study, mk.initial_M, concentrations[0]),
            cycles,
            [
                study.max_product_M(initial_concentrations(study, mk.initial_M, c))
                for c in concentrations
            ],
            list(mk.plot_species),
            mk.catalyst_snapshot_h,
            mk.conversion,
            study.output_step_s,
            product=study.product,
        )
        analysis.plot_catalyst_concentration_vs_reactant(log_x=True, figsize=(15, 10))
        analysis.plot_concentration_evolution(
            T, True, True, *grid(len(concentrations)), (20, 16)
        )
        analysis.plot_reactant_vs_product_conversion(log_x=True, figsize=(15, 10))
        shares = AuxiliaryFunctions.get_concentrations_of_catalyst(
            analysis.simulations_dfs, mk.catalyst_snapshot_h, cycles
        )
        thresholds = [m * mk.conversion for m in analysis.max_product_M]
        times = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
            analysis.simulations_dfs, concentrations, studied, thresholds, study.product
        )
        summary = {
            "concentrations_M": concentrations.tolist(),
            "catalyst_M": dict(zip(cycles, shares)),
            "t99_h": [t for t, _ in times],
        }
        Path("microkinetics.json").write_text(json.dumps(summary, indent=1) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="microkatc.py", description="Microkinetic analysis of catalytic cycles"
    )
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="run every analysis the study requests")
    run.add_argument("study")
    run.add_argument(
        "--fresh", action="store_true", help="delete the study's results/ first"
    )
    commands.add_parser(
        "check", help="check the study without running anything"
    ).add_argument("study")
    args = parser.parse_args(argv)

    try:
        study = load_study(args.study)
    except StudyError as error:
        print(f"{args.study}: invalid study", file=sys.stderr)
        for message in error.messages:
            print(f"  {message}", file=sys.stderr)
        return 1
    for warning in study.warnings:
        print(warning)
    if args.command == "check":
        print(
            f"{args.study}: OK ({len(study.steps)} steps, {len(study.species())} species)"
        )
        return 0

    try:
        results, info = prepare_results(study, args.fresh)
    except StaleResultsError as error:
        print(error, file=sys.stderr)
        return 1
    started, previous = time.time(), os.getcwd()
    os.chdir(results)
    try:
        run_analyses(study)
    except Exception as error:
        traceback.print_exc()
        print(f"{study.path}: run failed: {error}", file=sys.stderr)
        return 2
    finally:
        os.chdir(previous)
    info.update(
        finished=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        duration_s=round(time.time() - started, 1),
    )
    _write_run_info(results, info)
    print(f"Results in {results}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run the tests and the example**

```bash
python tests/test_microkatc.py
(cd /tmp && python "$OLDPWD/microkatc.py" run "$OLDPWD/examples/typed_energies/study.yaml"); echo "exit $?"
cat examples/typed_energies/results/microkinetics.json
rm -rf examples/typed_energies/results
```
Expected: `ok`; `exit 0` in about 3 s; three `t99_h` values that grow with the substrate concentration, and a catalyst total of 0.001 M in every column.

- [ ] **Step 6: Commit**

```bash
ruff check . && ruff format --check .
git add microkatc.py examples/typed_energies/study.yaml tests/test_microkatc.py
git commit -m "feat: microkatc.py run/check with results/, run_info.json and a stale-results guard"
```

---

### Task 8: Rewrite the hydroformylation example as a study file

**Files:**
- Create: `examples/hydroformylation/study.yaml`
- Move: `GaussOutputFiles/` → `examples/hydroformylation/GaussOutputFiles/`
- Delete: `reactions.csv`, `get_G_compounds.sh`, `bash_parsing.py`, `tests/test_g_step_failure.py` (its case is covered by `test_thermochemistry.py::test_zero_energy_from_thermochange_stops_the_run` and the file checks in `test_study.py`)
- Modify (patch below): `auxiliary_functions.py` (drop `calculate_G_values`, `find_intermediates_of_cycle`, `find_limiting_reactant_concentration`, `reactions_number`), `calculating_G_for_microkinetics.py` (keep `ReactionFileParser.parse_reaction`, `GibbsEnergyCalculator` and the constants only), `main.py` (one-line wrapper), `apparent_activation_energy.py` (messages no longer mention `reactions.csv`), `readme_figures.py` and `tests/check_reproduction.py` (read the example study and its `results/`; pass `time_step=STUDY.output_step_s` so cache names match the run), `tests/test_paper_barriers.py` (build the table from the study), `.gitignore` (`results/`), `.github/workflows/nightly.yml` (run the YAML example; run `check_equivalence.py`; run on pull requests that change code, the example or the checks)

**Interfaces:**
- Consumes: Tasks 4–7.
- Produces: `examples/hydroformylation/study.yaml`; `readme_figures.STUDY`, `MK`, `EA`, `REACTANT`, `T`, `CONVERSION`, `MK_SIMULATION_TIME`, `latest(kind)` (imported by `tests/check_reproduction.py`).

- [ ] **Step 1: Create the study file** (`examples/hydroformylation/study.yaml`)

The steps are `reactions.csv` row by row (`=` written as `<=>`, the TS column as `via`); the conditions are the current `main.py` values.

```yaml
# MicroKatc study: Rh-catalysed hydroformylation of ethylene with PMe3
# Abdullayev et al., ACS Catal. 2025, 15, 4739 (doi:10.1021/acscatal.5c00348)
name: Hydroformylation, 0L/1L ligand exchange

species:
  files: GaussOutputFiles          # <name>.out for every species and transition state
  product: prod
  overall_reaction: ete + CO + H2 <=> prod
  cycles:
    0L: [I1_0L, I2_0L, I3_0L, I4_0L, I5_0L, I6_0L, I7_0L, I8_0L, I9_0L]
    1L: [I1_1L, I2c_1L, I2t_1L, I3_1L, I4_1L, I5_1L, I6_1L, I7_1L, I8_1L, I9c_1L, I9t_1L]

steps:                             # every step is reversible; "via" names the transition state
  - I1_0L <=> I2_0L + CO
  - I2_0L + ete <=> I3_0L
  - I3_0L <=> I4_0L  via TS1_0L
  - I4_0L + CO <=> I5_0L
  - I5_0L <=> I6_0L  via TS2_0L
  - I6_0L + CO <=> I7_0L
  - I6_0L + H2 <=> I8_0L  via TS3_0L
  - I8_0L <=> I9_0L  via TS4_0L
  - I9_0L <=> prod + I2_0L
  - I2_0L + PMe3 <=> I1_1L
  - I1_1L <=> I2c_1L + CO
  - I2c_1L + ete <=> I3_1L
  - I3_1L <=> I4_1L  via TS1_1L
  - I4_1L + CO <=> I5_1L
  - I5_1L <=> I6_1L  via TS2_1L
  - I6_1L + CO <=> I7_1L
  - I6_1L + H2 <=> I8_1L  via TS3_1L
  - I8_1L <=> I9t_1L  via TS4t_1L
  - I8_1L <=> I9c_1L  via TS4c_1L
  - I9t_1L <=> I2t_1L + prod
  - I9c_1L <=> prod + I2c_1L
  - I2t_1L + CO <=> I1_1L

conditions:
  temperature_K: 350
  studied_species: PMe3
  studied_range_M: {from: 1e-10, to: 0.1, points: 19}

analyses:
  activation_energy:               # low catalyst concentration (differential reactor)
    initial_M: {CO: 0.05, H2: 0.05, ete: 0.05, I1_0L: 1e-6}
    temperatures_K: {from: 325, to: 375, points: 5}
    simulation_time_s: 10000
    sampling_time_h: 2
  degree_of_rate_control:          # uses the activation_energy initial concentrations
    e_shift_kcal_mol: 0.1
    simulation_time_s: 10000
    cores: 8
  microkinetics:
    initial_M: {CO: 0.05, H2: 0.05, ete: 0.05, I1_0L: 5e-4}
    simulation_time_s: 100000
    catalyst_snapshot_h: 1
    conversion: 0.99
    plot_species: [I7_0L, I7_1L, I1_0L, I1_1L]
```

- [ ] **Step 2: Move and delete the old inputs**

```bash
mkdir -p examples/hydroformylation
git mv GaussOutputFiles examples/hydroformylation/GaussOutputFiles
git rm reactions.csv get_G_compounds.sh bash_parsing.py tests/test_g_step_failure.py
python microkatc.py check examples/hydroformylation/study.yaml
```
Expected: `examples/hydroformylation/study.yaml: OK (22 steps, 25 species)` (`$thermochange` must be exported).

- [ ] **Step 3: Apply the patch**

Save the block below exactly as `task8.patch` in the repository root, then apply it:

```diff
--- a/auxiliary_functions.py
+++ b/auxiliary_functions.py
@@ -1,18 +1,7 @@
-"""Helper functions shared by the analyses (cycle intermediates, catalyst concentrations, G values, pressure)"""
+"""Helper functions shared by the analyses (cache names, catalyst concentrations, conversion times, pressure)"""
 
 import hashlib
 import json
-import os
-import re
-import subprocess
-
-import pandas as pd
-
-from calculating_G_for_microkinetics import (
-    G_COMPOUNDS_OUTPUT_DIR_NAME,
-    REACTION_DF_OUTPUT_DIR_NAME,
-)
-from file_operations import run_directory
 
 R_L_atm_per_mol_K = 0.082057366080960
 
@@ -29,34 +18,6 @@
             default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o),
         )
         return hashlib.md5(dumped.encode()).hexdigest()[:8]
-
-    @staticmethod
-    def find_intermediates_of_cycle(cycles):
-        """Returns intermediates belonging to a corresponding cycle as a dictionary"""
-        compounds = pd.read_csv("compounds.csv")["Compounds"].to_list()
-        cycles_intermediates_dict = {}
-
-        for cycle in cycles:
-            # Regex pattern for intermediates (I1_0L, I5_1L, I2c_1L, I9t_1L...)
-            pattern = re.compile(rf"I\d+[a-z]*_{cycle}")
-
-            # Find intermediates matching the pattern for the current cycle
-            intermediates = [
-                intermediate
-                for intermediate in compounds
-                if pattern.fullmatch(intermediate)
-            ]
-            cycles_intermediates_dict[cycle] = intermediates
-
-        return cycles_intermediates_dict
-
-    @staticmethod
-    def find_limiting_reactant_concentration(initial_concentrations, reactants):
-        """Returns concentration of a limiting reactant in a system"""
-        min_initial_reactant_concentrations = min(
-            [initial_concentrations[reactant] for reactant in reactants]
-        )
-        return min_initial_reactant_concentrations
 
     @staticmethod
     def get_concentration_for_cycle(simulations_dfs, time, intermediates):
@@ -123,53 +84,6 @@
         return times_conv_reac_conc
 
     @staticmethod
-    def reactions_number(reaction_df):
-        """Returns the reactions of the "Rx" column as a list"""
-        reactions = reaction_df["Rx"].to_list()
-        return reactions
-
-    @staticmethod
-    def calculate_G_values(temperature, pressure):
-        """Calculates G of compounds and reactions in a system at specified T and P"""
-        reaction_df_file_name = f"reaction_df_{temperature}K_{pressure:.5e}atm.csv"
-        G_compounds_file_name = f"G_values_at_{temperature}K_{pressure:.5e}atm.csv"
-
-        G_compounds_file_path = os.path.join(
-            run_directory(), G_COMPOUNDS_OUTPUT_DIR_NAME, G_compounds_file_name
-        )
-        reaction_df_file_path = os.path.join(
-            run_directory(), REACTION_DF_OUTPUT_DIR_NAME, reaction_df_file_name
-        )
-        outputs = [G_compounds_file_path, reaction_df_file_path]
-
-        # Calculations will be executed if either file for the specified T and P is missing
-        if not all(os.path.exists(path) for path in outputs):
-            # Path to the Bash script
-            script_path = os.path.join(run_directory(), "get_G_compounds.sh")
-            print(f"Creating and saving: {reaction_df_file_name}")
-            print(f"Creating and saving: {G_compounds_file_name}")
-            result = subprocess.run(
-                ["bash", script_path, f"{temperature}", f"{pressure}"],
-                capture_output=True,
-                text=True,
-                check=False,
-            )
-            # get_G_compounds.sh exits 0 even when thermochange or the Python step fails,
-            # so check for its output files instead of the exit code. A failure can leave the
-            # compound file without the reaction file; remove it so the next run starts clean.
-            missing = [path for path in outputs if not os.path.exists(path)]
-            if missing:
-                for path in outputs:
-                    if os.path.exists(path):
-                        os.remove(path)
-                raise RuntimeError(
-                    f"get_G_compounds.sh did not create {', '.join(missing)}. "
-                    "Is $thermochange exported, and does every species and transition state "
-                    "in reactions.csv have a .out file in GaussOutputFiles/?\n"
-                    f"{result.stdout}{result.stderr}"
-                )
-
-    @staticmethod
     def compute_pressure_value(temperature_value):
         """Computes pressure to satisfy liquid medium at standard conditions (1M)"""
         return 1 * R_L_atm_per_mol_K * temperature_value
--- a/calculating_G_for_microkinetics.py
+++ b/calculating_G_for_microkinetics.py
@@ -1,31 +1,12 @@
-"""Computes Gibbs energies of compounds and of direct/inverse reactions from thermochange output (called by get_G_compounds.sh)"""
-
-import os
-
-import pandas as pd
-
-from bash_parsing import BashParser
-from file_operations import FileOperations, run_directory
+"""Barriers of each step from the Gibbs energies of its species (used by thermochemistry.py)"""
 
 DIFFUSION_BARRIER = 4
 G_COMPOUNDS_OUTPUT_DIR_NAME = "G_values_of_compounds"
 REACTION_DF_OUTPUT_DIR_NAME = "G_values_of_reactions"
 
 
-class TemperaturePressureParser:
-    """Reads temperature and pressure from the command line"""
-
-    @staticmethod
-    def get_temperature_pressure():
-        """Fetch temperature and pressure values from bash input"""
-        return BashParser.get_temperature_pressure_from_input()
-
-
 class ReactionFileParser:
-    """Reads reactions.csv and splits reaction strings into reactants and products"""
-
-    def __init__(self, reaction_file="reactions.csv"):
-        self.reactions_df = pd.read_csv(reaction_file, sep=",")
+    """Splits reaction strings into reactants and products"""
 
     @staticmethod
     def parse_reaction(reaction):
@@ -98,87 +79,3 @@
             Ginv_list.append(Ginv)
 
         return Gdir_list, Ginv_list
-
-
-class CompoundsDataHandler:
-    """Writes compounds.csv and the compound G values for a given T and P"""
-
-    def __init__(self, df_compounds_file_name="compounds.csv"):
-        self.df_compounds_file_name = df_compounds_file_name
-
-    def create_compounds_df(self, compounds):
-        """Creates and saves dataframe with compounds in the system if it doesn't exist"""
-        compounds_df_path = os.path.join(run_directory(), self.df_compounds_file_name)
-        if not os.path.exists(compounds_df_path):
-            print("Creating compounds.csv file...")
-            compounds.to_csv(f"{self.df_compounds_file_name}", index=False)
-
-    def save_compounds_G_values(self, compound_energy_dataframe, G_compounds_file_name):
-        """Saves compounds data to CSV and move to the appropriate directory"""
-        os.makedirs(G_COMPOUNDS_OUTPUT_DIR_NAME, exist_ok=True)
-        compound_energy_dataframe.to_csv(G_compounds_file_name, index=True)
-        FileOperations.move_to_output_directory(
-            G_COMPOUNDS_OUTPUT_DIR_NAME, G_compounds_file_name
-        )
-
-
-class ReactionDataHandler:
-    """Writes the reaction dataframe with Gdir and Ginv for a given T and P"""
-
-    def save_reaction_data(self, reactions_df, reaction_df_file_name):
-        """Saves reaction data to CSV and move to the appropriate directory"""
-        os.makedirs(REACTION_DF_OUTPUT_DIR_NAME, exist_ok=True)
-        reactions_df.to_csv(reaction_df_file_name, index=False)
-        FileOperations.move_to_output_directory(
-            REACTION_DF_OUTPUT_DIR_NAME, reaction_df_file_name
-        )
-
-
-class ReactionCalculator:
-    """Runs the G calculation for the T and P given on the command line and saves the results"""
-
-    def __init__(self, reaction_file="reactions.csv"):
-        # Initialize parsers and handlers
-        self.temperature, self.pressure = (
-            TemperaturePressureParser.get_temperature_pressure()
-        )
-        self.reaction_file_parser = ReactionFileParser(reaction_file)
-        self.compounds_data_handler = CompoundsDataHandler()
-        self.reaction_data_handler = ReactionDataHandler()
-
-        self.G_compounds_file_name = (
-            f"G_values_at_{self.temperature}K_{self.pressure:.5e}atm.csv"
-        )
-        self.reaction_df_file_name = (
-            f"reaction_df_{self.temperature}K_{self.pressure:.5e}atm.csv"
-        )
-
-        # Fetch compound data and calculate Gibbs free energies
-        self.compound_energy_dataframe, self.compounds = (
-            BashParser.parse_output_from_bash_script()
-        )
-        self.compounds_data_handler.create_compounds_df(self.compounds)
-        self.compounds_data_handler.save_compounds_G_values(
-            self.compound_energy_dataframe, self.G_compounds_file_name
-        )
-
-        # Calculate Gibbs free energy for reactions
-        self.gibbs_calculator = GibbsEnergyCalculator(self.compound_energy_dataframe)
-        Gdir_values, Ginv_values = (
-            self.gibbs_calculator.calculate_direct_inverse_reactions_gibbs_free_energies(
-                self.reaction_file_parser.reactions_df
-            )
-        )
-
-        # Update reaction dataframe and save
-        (
-            self.reaction_file_parser.reactions_df["Gdir"],
-            self.reaction_file_parser.reactions_df["Ginv"],
-        ) = Gdir_values, Ginv_values
-        self.reaction_data_handler.save_reaction_data(
-            self.reaction_file_parser.reactions_df, self.reaction_df_file_name
-        )
-
-
-if __name__ == "__main__":
-    ReactionCalculator()
--- a/main.py
+++ b/main.py
@@ -1,207 +1,15 @@
-"""Entry point: apparent Ea, DRC and concentration evolution analyses of the cycle in reactions.csv"""
+"""Runs the hydroformylation example of the paper (ACS Catal. 2025, 15, 4739).
 
-import matplotlib.pyplot as plt
-import numpy as np
-import pandas as pd
+For your own system, write a study file and run: python microkatc.py run <study.yaml>
+"""
 
-from apparent_activation_energy import ApparentEaAnalysis, DRCAnalysis
-from auxiliary_functions import AuxiliaryFunctions
-from microkinetics_simulation import MicroKinetics
+import sys
+from pathlib import Path
 
-
-def main():
-    """Runs the apparent Ea, DRC and microkinetics analyses with the parameters set below"""
-    # Load reaction data
-    reaction_df = pd.read_csv("reactions.csv", sep=",")
-    reactions = AuxiliaryFunctions.reactions_number(reaction_df)
-
-    # Specify temperature in K
-    temperature_value = 350.0
-
-    # Specify range of T for which simulations will be calculated (Ea analysis)
-    T_values_array_Ea = np.linspace(
-        start=temperature_value - 25, stop=temperature_value + 25, num=5, endpoint=True
-    )
-
-    # Concentrations of the studied reactant (log scale): every half decade from 1e-10 to 0.1 M, as in the paper
-    left_border_concentration = -10
-    right_border_concentration = -1
-    reactant_concentration_array = np.logspace(
-        start=left_border_concentration,
-        stop=right_border_concentration,
-        num=19,
-        endpoint=True,
-    )
-
-    # Initial concentrations (M) for the apparent Ea and DRC analyses: a very low catalyst
-    # concentration, as in the paper (ACS Catal. 2025, 15, 4739)
-    c0_1 = {
-        "CO": 0.05,
-        "H2": 0.05,
-        "ete": 0.05,
-        "prod": 0,
-        "I1_0L": 1e-6,
-        "PMe3": 0.0000005,
-    }
-
-    # Specify a reactant to study
-    reactant_to_study = "PMe3"
-
-    # Specify simulation time for apparent Ea analysis (in seconds)
-    total_simulation_time_Ea = 10_000
-
-    # Specify time in h at which rate of reactions (fluxes) is considered
-    time = 2
-
-    # Specify time step for the microkinetics simulations
-    time_step = 1
-
-    # Barrier tables for every temperature, before any simulation (the analyses no longer compute them)
-    for T_value in sorted({*T_values_array_Ea, temperature_value}):
-        AuxiliaryFunctions.calculate_G_values(
-            T_value, AuxiliaryFunctions.compute_pressure_value(T_value)
-        )
-
-    # Perform MicroKinetics Analysis
-    analysis1 = ApparentEaAnalysis(
-        temperature_value,
-        T_values_array_Ea,
-        reactant_concentration_array,
-        reactant_to_study,
-        c0_1,
-        total_simulation_time_Ea,
-        time,
-        reactions,
-        time_step,
-    )
-
-    # Specify reactant initial concentration from the concentration array that will be considered for ln(ri) vs. 1/T plot
-    reactant_initial_concentration_plot = reactant_concentration_array[0]
-
-    # Specify reactions that will be on plot of Ea vs. c0 (flux based)
-    reactions_to_plot = [reaction for reaction in reactions if "prod" in reaction]
-
-    # Specify compounds that will be on plot of Ea vs. c0 (compound based)
-    compounds_to_plot = ["prod"]
-
-    # Calculates and saves as csv dataframe with ri (flux) related parameters if wasn't calculated before for the specified parameters
-    df_flux = analysis1.df_flux
-
-    # Calculates and saves as csv dataframe with vi (rate) related parameters if wasn't calculated before for the specified parameters
-    df_rate = analysis1.df_rate
-
-    # Adjust nrows and ncols to your case
-    analysis1.plot_ln_ri_vs_1_over_T(
-        df_flux=df_flux,
-        reactant_initial_concentration_plot=reactant_initial_concentration_plot,
-        nrows=5,
-        ncols=5,
-        figsize=(15, 15),
-    )
-    analysis1.plot_Ea_vs_c0_flux_based(
-        df_flux=df_flux,
-        reactions_to_plot=reactions_to_plot,
-        nrows=2,
-        ncols=3,
-        figsize=(15, 10),
-        log_x=True,
-    )
-    analysis1.plot_Ea_vs_c0_rate_based(
-        df_rate=df_rate,
-        compounds_to_plot=compounds_to_plot,
-        nrows=1,
-        ncols=1,
-        figsize=(15, 10),
-        log_x=True,
-    )
-
-    # Barrier shift (kcal/mol) for the central-difference degree of rate control (DRC). The paper
-    # used a one-sided 0.01 shift; +-0.1 averaged is 8 times less noisy and agrees with it within 0.035
-    e_shift = 0.1
-
-    # Temperature(s) at which DRC is calculated: the working temperature, as in the paper
-    T_values_array_drc = np.array([temperature_value])
-
-    # Specify simulation time for apparent Ea analysis (in seconds)
-    total_simulation_time_drc = 10_000
-
-    # Specify number of cores to compute drc
-    cores = 8
-
-    analysis2 = DRCAnalysis(
-        reactant_concentration_array,
-        T_values_array_drc,
-        total_simulation_time_drc,
-        c0_1,
-        reactant_to_study,
-        reactions,
-        e_shift,
-        cores,
-    )
-
-    # Temperature(s) that will be considered for drc vs. c0 plot (this variable can be an array if multiple plots are needed)
-    temperature_value_drc_plot = T_values_array_drc[0]
-
-    # Calculates and saves as csv dataframe with degree rate control (drc) related parameters if wasn't calculated before for the specified parameters
-    df_drc = analysis2.df_drc
-
-    # Adjust nrows and ncols to your case
-    analysis2.plot_drc_vs_c0(
-        nrows=5,
-        ncols=5,
-        figsize=(15, 15),
-        df_drc=df_drc,
-        temperature_value=temperature_value_drc_plot,
-        log_x=True,
-    )
-
-    # Initial concentrations (M) for the catalyst distribution and conversion time, as in the paper
-    c0_2 = {"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 5e-4, "PMe3": 1}
-
-    # Specify cycle names
-    cycles = ["0L", "1L"]
-
-    # Specify simulation time for "normal" microkinetics (in seconds)
-    total_simulation_time_microkinetics = 100_000
-
-    # Define reactants that could affect final concentration of a product (can be limiting reactant)
-    main_reactants = ["CO", "H2", "ete"]
-
-    # Specify compounds that will appear on a concentration evolution plot
-    compounds_to_plot = ["I7_0L", "I7_1L", "I1_0L", "I1_1L"]
-
-    # Set the time (in h) at which concentration of catalyst will be calculated (for C(catalyst) in each cycle in a system vs. initial concentration of studied reactant)
-    catalyst_concentration_time = 1
-
-    # Set the percentage of product conversion
-    percentage_of_convertion = 0.99
-
-    # Temperature of the concentration evolution plot: the simulations above run at temperature_value
-    temperatures = temperature_value
-
-    analysis3 = MicroKinetics(
-        temperature_value,
-        total_simulation_time_microkinetics,
-        reactant_concentration_array,
-        reactant_to_study,
-        c0_2,
-        AuxiliaryFunctions.find_intermediates_of_cycle(cycles),
-        [min(c0_2[r] for r in main_reactants)] * len(reactant_concentration_array),
-        compounds_to_plot,
-        catalyst_concentration_time,
-        percentage_of_convertion,
-        time_step,
-    )
-
-    # Adjust nrows and ncols to your case
-    analysis3.plot_catalyst_concentration_vs_reactant(log_x=True, figsize=(15, 10))
-    analysis3.plot_concentration_evolution(
-        temperatures, log_y=True, log_x=True, nrows=4, ncols=5, figsize=(20, 16)
-    )
-    analysis3.plot_reactant_vs_product_conversion(log_x=True, figsize=(15, 10))
-
-    plt.show()
-
+from microkatc import main
 
 if __name__ == "__main__":
-    main()
+    example = (
+        Path(__file__).resolve().parent / "examples" / "hydroformylation" / "study.yaml"
+    )
+    sys.exit(main(["run", str(example), *sys.argv[1:]]))
--- a/apparent_activation_energy.py
+++ b/apparent_activation_energy.py
@@ -116,11 +116,11 @@
 
 
 def reaction_of_flux(flux_name, reactions):
-    """Returns the reactions.csv row of a COPASI flux column: copasi_helper names row i r{i+1:02d}"""
+    """Returns the step of a COPASI flux column: copasi_helper names step i r{i+1:02d}"""
     match = re.fullmatch(r"r(\d+)\.Flux", flux_name)
     if match is None or not 1 <= int(match.group(1)) <= len(reactions):
         raise ValueError(
-            f"Cannot match flux column {flux_name!r} to a row of reactions.csv"
+            f"Cannot match flux column {flux_name!r} to a step of the study"
         )
     return reactions[int(match.group(1)) - 1]
 
@@ -145,10 +145,10 @@
             else [key for key in df["name"].unique() if key.endswith(".Rate")]
         )
 
-        # Every reactions.csv row must have exactly one flux column
+        # Every step must have exactly one flux column
         if calculation_type == "ri" and len(keys) != len(reactions):
             raise ValueError(
-                f"COPASI returned {len(keys)} '.Flux' columns but reactions.csv has "
+                f"COPASI returned {len(keys)} '.Flux' columns but the study has "
                 f"{len(reactions)} reactions; cannot match fluxes to reactions"
             )
 
@@ -488,11 +488,11 @@
                     axis=0,
                 )
 
-                # DRC coefficients are matched to reactions.csv rows by position
+                # DRC coefficients are matched to the steps by position (drc_calc follows step order)
                 if len(drc_coefficients) != len(self.reactions):
                     raise ValueError(
                         f"COPASI returned {len(drc_coefficients)} DRC coefficients but "
-                        f"reactions.csv has {len(self.reactions)} reactions; cannot match them"
+                        f"the study has {len(self.reactions)} steps; cannot match them"
                     )
 
                 data.append(
--- a/readme_figures.py
+++ b/readme_figures.py
@@ -1,31 +1,40 @@
-"""Builds the README figures in pics/ from the results that main.py saves.
+"""Builds the README figures in pics/ from the results of the hydroformylation example.
 
 Each figure is sized for GitHub's README column (about 840 px wide) and shares one
 style: colour identifies the catalytic cycle (0L blue, 1L pink, as in the paper's TOC
 graphic) and the product is green with a dashed line and triangles. main.py's own
 figures, with every step and every concentration, stay in microkinetics_simulations_images/.
 
-Run after main.py, from the same directory:
+Run after the example (python main.py, or python microkatc.py run examples/hydroformylation/study.yaml):
     python readme_figures.py
 """
 
 import glob
 import os
+from pathlib import Path
 
 import matplotlib.pyplot as plt
 import numpy as np
 import pandas as pd
 
 from auxiliary_functions import AuxiliaryFunctions
+from microkatc import initial_concentrations
 from microkinetics_simulation import SIMULATIONS_OUTPUT_DIR_NAME, SimulationHandler
-
-REACTANT = "PMe3"
-T = 350.0
-# Same as c0_1 / c0_2 and the simulation times in main.py, so the cached simulations are reused
-C0_EA = {"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 1e-6}
-C0_MK = {"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 5e-4}
-EA_SIMULATION_TIME, EA_TIME_H = 10_000, 2
-MK_SIMULATION_TIME, CATALYST_TIME_H, CONVERSION = 100_000, 1, 0.99
+from study import load_study
+
+REPO = Path(__file__).resolve().parent
+# Every condition comes from the example study, so the cached simulations of its run are reused
+STUDY = load_study(
+    REPO / "examples" / "hydroformylation" / "study.yaml", check_energy_sources=False
+)
+REACTANT, T, PRODUCT_NAME = STUDY.studied_species, STUDY.temperature_K, STUDY.product
+EA, MK = STUDY.activation_energy, STUDY.microkinetics
+EA_SIMULATION_TIME, EA_TIME_H = EA.simulation_time_s, EA.sampling_time_h
+MK_SIMULATION_TIME, CATALYST_TIME_H, CONVERSION = (
+    MK.simulation_time_s,
+    MK.catalyst_snapshot_h,
+    MK.conversion,
+)
 
 CYCLE_0L, CYCLE_1L, PRODUCT, INK, MUTED = (
     "#2a78d6",
@@ -88,7 +97,7 @@
 
 
 def save(fig, name):
-    fig.savefig(os.path.join("pics", name), dpi=170, bbox_inches="tight")
+    fig.savefig(REPO / "pics" / name, dpi=170, bbox_inches="tight")
     plt.close(fig)
     print(f"Saved pics/{name}")
 
@@ -99,6 +108,7 @@
 
 
 def main():
+    os.chdir(STUDY.results_dir)
     flux, rate, drc = latest("flux"), latest("rate"), latest("drc")
     flux["reaction"] = flux["reaction"].str.split().str.join(" ")
     drc = drc.rename(columns=lambda c: " ".join(c.split())).sort_values(
@@ -188,10 +198,16 @@
 
     # 4-6. Catalyst distribution, concentration profiles and conversion time (main.py's c0_2)
     mk = SimulationHandler(
-        T, AuxiliaryFunctions.compute_pressure_value(T), REACTANT, MK_SIMULATION_TIME
-    )
-    sims = [mk.get_simulation_df({**C0_MK, REACTANT: c}) for c in c0]
-    cycles = AuxiliaryFunctions.find_intermediates_of_cycle(["0L", "1L"])
+        T,
+        AuxiliaryFunctions.compute_pressure_value(T),
+        REACTANT,
+        MK_SIMULATION_TIME,
+        STUDY.output_step_s,
+    )
+    sims = [
+        mk.get_simulation_df(initial_concentrations(STUDY, MK.initial_M, c)) for c in c0
+    ]
+    cycles = {label: list(members) for label, members in STUDY.cycles.items()}
     catalyst = AuxiliaryFunctions.get_concentrations_of_catalyst(
         sims, CATALYST_TIME_H, cycles
     )
@@ -202,7 +218,7 @@
         catalyst, (CYCLE_0L, CYCLE_1L), ("0L cycle", "1L cycle")
     ):
         ax.plot(c0, 100 * np.array(share) / total, "o-", color=colour, label=label)
-    ax.set_ylabel(f"share of catalyst at t = {CATALYST_TIME_H} h (%)")
+    ax.set_ylabel(f"share of catalyst at t = {CATALYST_TIME_H:g} h (%)")
     ax.set_ylim(-3, 103)
     ax.set_title("Catalyst distribution between the cycles", color=INK)
     log_x(ax)
@@ -235,7 +251,15 @@
     save(fig, "concentration_profiles.png")
 
     times = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
-        sims, c0, REACTANT, CONVERSION * min(C0_MK[r] for r in ("CO", "H2", "ete"))
+        sims,
+        c0,
+        REACTANT,
+        [
+            CONVERSION
+            * STUDY.max_product_M(initial_concentrations(STUDY, MK.initial_M, c))
+            for c in c0
+        ],
+        PRODUCT_NAME,
     )
     fig, ax = plt.subplots(figsize=(8.5, 4.3))
     ax.plot(
@@ -253,10 +277,16 @@
 
     # 7. Poisoning intermediates at the Ea sampling time, against [PMe3] (main.py's c0_1)
     ea = SimulationHandler(
-        T, AuxiliaryFunctions.compute_pressure_value(T), REACTANT, EA_SIMULATION_TIME
+        T,
+        AuxiliaryFunctions.compute_pressure_value(T),
+        REACTANT,
+        EA_SIMULATION_TIME,
+        STUDY.output_step_s,
     )
     rows = [
-        ea.get_simulation_df({**C0_EA, REACTANT: c}).set_index("time").loc[EA_TIME_H]
+        ea.get_simulation_df(initial_concentrations(STUDY, EA.initial_M, c))
+        .set_index("time")
+        .loc[EA_TIME_H]
         for c in c0
     ]
     fig, ax = plt.subplots(figsize=(8.5, 4.6))
@@ -274,7 +304,7 @@
             label=species.replace("_", "-"),
         )
     ax.set(yscale="log", ylabel="concentration (M)")
-    ax.set_title(f"Poisoning intermediates at t = {EA_TIME_H} h", color=INK)
+    ax.set_title(f"Poisoning intermediates at t = {EA_TIME_H:g} h", color=INK)
     log_x(ax)
     ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
     save(fig, "poisoning_intermediates.png")
--- a/tests/check_reproduction.py
+++ b/tests/check_reproduction.py
@@ -1,7 +1,6 @@
 """Checks a full main.py run against the published results (Abdullayev et al., ACS Catal. 2025, 15, 4739).
 
-Not a unit test: it reads the tables and simulations that main.py saves, so run it after main.py,
-from the repository root:
+Not a unit test: it reads the results of the example study, so run it after the example:
     python main.py && python tests/check_reproduction.py
 
 Paper values come from the text and from reading Figures 4-6; tolerances cover that reading
@@ -16,15 +15,19 @@
 sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
 
 from auxiliary_functions import AuxiliaryFunctions
+from microkatc import initial_concentrations
 from microkinetics_simulation import SimulationHandler
 from readme_figures import (
-    C0_MK,
     CONVERSION,
+    MK,
     MK_SIMULATION_TIME,
     REACTANT,
+    STUDY,
     T,
     latest,
 )
+
+os.chdir(STUDY.results_dir)
 
 failures = []
 
@@ -86,10 +89,17 @@
 
 # Figure 4: catalyst distribution at 1 h and time to 99 % conversion
 handler = SimulationHandler(
-    T, AuxiliaryFunctions.compute_pressure_value(T), REACTANT, MK_SIMULATION_TIME
+    T,
+    AuxiliaryFunctions.compute_pressure_value(T),
+    REACTANT,
+    MK_SIMULATION_TIME,
+    STUDY.output_step_s,
 )
-sims = [handler.get_simulation_df({**C0_MK, REACTANT: c}) for c in c0]
-cycles = AuxiliaryFunctions.find_intermediates_of_cycle(["0L", "1L"])
+sims = [
+    handler.get_simulation_df(initial_concentrations(STUDY, MK.initial_M, c))
+    for c in c0
+]
+cycles = {label: list(members) for label, members in STUDY.cycles.items()}
 zero_l, one_l = AuxiliaryFunctions.get_concentrations_of_catalyst(sims, 1, cycles)
 crossing = c0[np.argmax(np.array(one_l) > np.array(zero_l))]
 check(
@@ -100,7 +110,15 @@
 )
 times = (
     AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
-        sims, c0, REACTANT, CONVERSION * min(C0_MK[r] for r in ("CO", "H2", "ete"))
+        sims,
+        c0,
+        REACTANT,
+        [
+            CONVERSION
+            * STUDY.max_product_M(initial_concentrations(STUDY, MK.initial_M, c))
+            for c in c0
+        ],
+        STUDY.product,
     )
 )
 check("time to 99 % conversion, low PMe3 / h (paper 5.66)", times[0][0], 5.55, 5.75)
--- a/tests/test_paper_barriers.py
+++ b/tests/test_paper_barriers.py
@@ -1,19 +1,16 @@
 """Gibbs barriers at 350 K must match Table S3 of the paper's Supporting Information.
 
-Abdullayev et al., ACS Catal. 2025, 15, 4739 (doi:10.1021/acscatal.5c00348). Runs the real
-thermochange step on GaussOutputFiles/, so it needs $thermochange exported; it is skipped otherwise.
+Abdullayev et al., ACS Catal. 2025, 15, 4739 (doi:10.1021/acscatal.5c00348). Builds the barrier
+table of examples/hydroformylation/study.yaml with the real thermochange, so it needs $thermochange
+exported; it is skipped otherwise.
     python tests/test_paper_barriers.py   (or: pytest tests)
 """
 
 import os
-import shutil
-import subprocess
 import sys
-import tempfile
-
-import pandas as pd
 
 ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
+sys.path.insert(0, ROOT)
 
 # Rx, Gdir, Ginv in kcal/mol at 350.0 K and 1.0 M (SI Table S3, rounded to 0.1)
 TABLE_S3 = [
@@ -46,18 +43,11 @@
     if not os.environ.get("thermochange"):
         print("skipped: export thermochange=/path/to/thermochange to run")
         return
-    work = os.path.join(tempfile.mkdtemp(), "MicroKatc")
-    shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns(".git", "tests", "pics"))
-    code = (
-        "from auxiliary_functions import AuxiliaryFunctions as A;"
-        "A.calculate_G_values(350.0, A.compute_pressure_value(350.0))"
-    )
-    subprocess.run([sys.executable, "-c", code], cwd=work, check=True)
+    from study import load_study
+    from thermochemistry import barrier_table, gibbs_energies
 
-    out = os.path.join(
-        work, "G_values_of_reactions", "reaction_df_350.0K_2.87201e+01atm.csv"
-    )
-    df = pd.read_csv(out)
+    study = load_study(os.path.join(ROOT, "examples", "hydroformylation", "study.yaml"))
+    df = barrier_table(study, gibbs_energies(study, 350.0))
     assert len(df) == len(TABLE_S3)
     for (rx, gdir, ginv), (_, row) in zip(TABLE_S3, df.iterrows()):
         assert " ".join(row["Rx"].split()) == rx
--- a/.gitignore
+++ b/.gitignore
@@ -160,3 +160,4 @@
 #  and can be added to the global gitignore or merged into this file.  For a more nuclear
 #  option (not recommended) you can uncomment the following to ignore the entire idea folder.
 #.idea/
+results/
--- a/.github/workflows/nightly.yml
+++ b/.github/workflows/nightly.yml
@@ -7,11 +7,16 @@
   schedule:
     - cron: "17 3 * * *"
   workflow_dispatch:
-  # Also when a pull request changes this workflow or the check, so edits are tested before merging
+  # Also when a pull request changes the code, the example or the checks, so the full results are
+  # compared with the paper and with the pre-YAML baseline before merging
   pull_request:
     paths:
+      - "*.py"
+      - requirements.txt
+      - examples/hydroformylation/**
+      - tests/check_*.py
+      - tests/data/**
       - .github/workflows/nightly.yml
-      - tests/check_reproduction.py
 
 permissions:
   contents: read
@@ -46,7 +51,7 @@
       - name: Run the full pipeline
         env:
           MPLBACKEND: Agg
-        run: python main.py
+        run: python microkatc.py run examples/hydroformylation/study.yaml
 
       - name: Build the README figures
         env:
@@ -58,6 +63,11 @@
           set -o pipefail
           python tests/check_reproduction.py | tee -a "$GITHUB_STEP_SUMMARY"
 
+      - name: Check results against the pre-YAML baseline
+        run: |
+          set -o pipefail
+          python tests/check_equivalence.py examples/hydroformylation/results | tee -a "$GITHUB_STEP_SUMMARY"
+
       - name: Keep figures and result tables
         if: always()
         uses: actions/upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a # v7.0.1
@@ -67,6 +77,8 @@
           if-no-files-found: warn
           path: |
             pics/*.png
-            microkinetics_simulations_images/
-            microkinetics_simulations/df_*.csv
-            G_values_of_reactions/
+            examples/hydroformylation/results/run_info.json
+            examples/hydroformylation/results/microkinetics.json
+            examples/hydroformylation/results/microkinetics_simulations_images/
+            examples/hydroformylation/results/microkinetics_simulations/df_*.csv
+            examples/hydroformylation/results/G_values_of_reactions/
```

```bash
git apply --check task8.patch && git apply task8.patch && rm task8.patch
grep -rn "find_intermediates_of_cycle\|calculate_G_values\|get_G_compounds\|bash_parsing" --include=*.py . || echo "old pipeline gone"
```
Expected: `old pipeline gone`.

- [ ] **Step 4: Run every test, including the SI Table S3 check**

```bash
for t in tests/test_*.py; do echo "$t: $(python $t 2>&1 | tail -1)"; done
```
Expected: `ok` for every file. `test_paper_barriers.py` must print `ok`, not `skipped`, so export `thermochange` first.

- [ ] **Step 5: Commit**

```bash
ruff check . && ruff format --check .
git add -A
git commit -m "feat: the hydroformylation example as examples/hydroformylation/study.yaml"
```

---

### Task 9: Stoichiometry against the analytical solution

**Files:**
- Create: `tests/test_stoichiometry.py`

**Interfaces:**
- Consumes: `copasi_parser` (copasi_helper `b9f26fe4`); needs COPASI, which CI installs.

- [ ] **Step 1: Write the test** (`tests/test_stoichiometry.py`)

```python
"""COPASI follows the rate law and mass balance of a step with a coefficient (real COPASI).

For 2 A -> B (reverse barrier so high that the step is irreversible), mass action gives
rate = k[A]^2 with two A consumed per event, so A(t) = A0 / (1 + 2 k A0 t) and A + 2 B = A0.
Both spellings MicroKatc can hand to copasi_helper ("A + A = B" and "2*A = B") must follow it.
    python tests/test_stoichiometry.py
"""

import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

T, A0, GDIR, GINV, END = 350.0, 0.1, 20.0, 60.0, 2000


def simulate(equation):
    import copasi_parser as cpx

    os.chdir(tempfile.mkdtemp())
    pd.DataFrame(
        {"Rx": [equation], "TS": ["-"], "Gdir": [GDIR], "Ginv": [GINV]}
    ).to_csv("rx.csv", index=False)
    model = cpx.prepare_copasi_model(
        reactions_file="rx.csv",
        temp=T,
        initial_concentrations={"A": A0, "B": 0},
        csv_delim=",",
    )
    cpx.time_course_simulation(model, END, 1, outfile="sim.txt")
    model.self_destruct()
    return cpx.read_simulation("sim.txt"), cpx.k_calc_Eyring(GDIR, T)


def test_second_order_rate_law_and_mass_balance():
    for equation in ("A + A = B", "2*A = B"):
        df, k = simulate(equation)
        t, a, b = df["time"].to_numpy(), df["A"].to_numpy(), df["B"].to_numpy()
        analytic = A0 / (1 + 2 * k * A0 * t)
        assert np.abs(a - analytic).max() / A0 < 1e-5, (
            equation,
            np.abs(a - analytic).max() / A0,
        )
        # COPASI writes 6 significant digits, so A + 2 B rounds to within 1e-7 M (measured 1.0e-7)
        assert np.abs(a + 2 * b - A0).max() < 5e-7, equation
        first_order = A0 * np.exp(
            -k * t
        )  # the law a missing coefficient would give: must differ
        assert np.abs(a - first_order).max() / A0 > 0.1, equation


if __name__ == "__main__":
    test_second_order_rate_law_and_mass_balance()
    print("ok")
```

- [ ] **Step 2: Run it**

Run: `python tests/test_stoichiometry.py`
Expected: `ok`. Measured on 2026-10-02: |A − analytic| / A0 = 8.3 × 10<sup>-7</sup>; |A + 2B − A0| = 1.0 × 10<sup>-7</sup> M, set by COPASI's 6 significant output digits.

- [ ] **Step 3: Negative control**

```bash
sed 's/A + A = B", "2\*A = B"/A = B", "A = B"/' tests/test_stoichiometry.py > /tmp/wrong.py && python /tmp/wrong.py; echo "exit $?"
```
Expected: an `AssertionError` and `exit 1` (a first-order step must not pass).

- [ ] **Step 4: Commit**

```bash
ruff check . && ruff format --check .
git add tests/test_stoichiometry.py
git commit -m "test: COPASI follows second-order kinetics and mass balance for 2 A <=> B"
```

---

### Task 10: Documentation

**Files:**
- Replace: `docs/USAGE.md` (content below)
- Modify: `README.md`, `CLAUDE.md` (patch below)

- [ ] **Step 1: Replace `docs/USAGE.md`**

```markdown
# Using MicroKatc on your own system

You describe a system in one **study file** (YAML) and run it with one command. Two complete examples ship with the repository:

- [`examples/hydroformylation/study.yaml`](../examples/hydroformylation/study.yaml): the paper's system, with energies from Gaussian output files. Run it once before writing your own study, to confirm your installation reproduces the published numbers.
- [`examples/typed_energies/study.yaml`](../examples/typed_energies/study.yaml): a small catalytic cycle with typed Gibbs energies, for when you have energies but no Gaussian files. It runs in a few seconds.

Contents:
- [1. Install](#1-install)
- [2. Write a study file](#2-write-a-study-file)
- [3. Check and run](#3-check-and-run)
- [4. Outputs](#4-outputs)
- [5. Troubleshooting](#5-troubleshooting)

## 1. Install

You need Python 3.10 and Bash (Linux or macOS; Windows through WSL). thermochange is only needed for studies that read output files.

```bash
git clone https://gitlab.com/dgarayr/thermochange.git
git clone https://github.com/0rkhann/MicroKatc.git
cd MicroKatc
python3.10 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt          # includes COPASI's Python bindings and copasi_helper
export thermochange=/path/to/thermochange
```

`thermochange` must be **exported**, not just set: MicroKatc runs it in subprocesses, which only see exported variables.

Check the installation against the paper:

```bash
python microkatc.py run examples/hydroformylation/study.yaml
python tests/check_reproduction.py
```

The last line should read `All results match the paper.`

## 2. Write a study file

A study file has four parts: `species`, `steps`, `conditions` and `analyses`. Start from one of the examples.

### 2.1 Steps

One line per elementary step. Every step is reversible.

```yaml
steps:
  - I1_0L <=> I2_0L + CO           # barrierless
  - I3_0L <=> I4_0L  via TS1_0L    # transition state named after "via"
  - 2 A <=> A2  via TS2            # coefficients: "2 A" or "2*A"
```

- `<=>`, `=` and `⇌` all mean the same.
- Species are separated by ` + `. Names are case-sensitive, start with a letter and contain no spaces.
- A **barrierless** step (no `via`) is treated as diffusion-controlled: its transition state is placed 4 kcal mol<sup>-1</sup> above the higher of its two sides (Besora et al., 2018).
- A **coefficient** n multiplies the species' Gibbs energy in the barrier, G(TS) − n·G(A), and gives a rate law of order n in that species.

### 2.2 Species and energies

```yaml
species:
  files: GaussOutputFiles          # either: one <name>.out per species and transition state
  product: prod
  overall_reaction: ete + CO + H2 <=> prod
  cycles:
    0L: [I1_0L, I2_0L, ...]
    1L: [I1_1L, I2c_1L, ...]
```

| Key | Meaning |
| --- | --- |
| `files` | Folder, relative to the study file, with the Gaussian (or ADF) output of a frequency calculation for every species and transition state. The file name without `.out` is the species name. |
| `energies_kcal_mol` | Instead of `files`: typed Gibbs energies (section 2.3). A study uses one or the other. |
| `product` | The species whose formation the degree of rate control and the conversion time follow. |
| `overall_reaction` | The net reaction, with coefficients. It sets the maximum yield: the smallest initial concentration ÷ coefficient over its reactants, times the product's coefficient. |
| `cycles` | Cycle label → the catalyst intermediates in that cycle. Needed for the `microkinetics` analysis. |
| `vibrational_correction` | `RRHO` (default) or `Grimme` (quasi-RRHO, thermochange's `-g`). Output files only. |

Output-file energies are corrected by thermochange to each temperature and a 1 M standard state.

### 2.3 Typed energies

```yaml
species:
  energies_kcal_mol:
    I1_0L: 0.0                                     # used at the working temperature only
    TS1_0L: {325: 11.5, 337.5: 11.5, 350: 11.5, 362.5: 11.5, 375: 11.5}
```

- Gibbs energies in kcal mol<sup>-1</sup> at 1 M, on any common reference: only differences enter the barriers.
- A single number is used at the working temperature only. The `activation_energy` analysis needs a value at every one of its temperatures, because the Arrhenius fit depends on how G changes with temperature.

### 2.4 Conditions

```yaml
conditions:
  temperature_K: 350
  studied_species: PMe3
  studied_range_M: {from: 1e-10, to: 0.1, points: 19}   # log-spaced
  output_step_s: 1                                      # optional, default 1
```

Every analysis is repeated for each concentration of the studied species. `output_step_s` is the time between output points of every simulation; sampling and snapshot times must fall on that grid.

### 2.5 Analyses

Request any of the three; those left out are skipped.

```yaml
analyses:
  activation_energy:               # needs a low catalyst concentration (differential reactor)
    initial_M: {CO: 0.05, H2: 0.05, ete: 0.05, I1_0L: 1e-6}
    temperatures_K: {from: 325, to: 375, points: 5}       # at least 3
    simulation_time_s: 10000
    sampling_time_h: 2             # rates are read at this time for the Arrhenius fits
  degree_of_rate_control:
    e_shift_kcal_mol: 0.1          # central difference: every barrier shifted up and down
    simulation_time_s: 10000       # rates are read at half of it
    cores: 8                       # optional, default: all CPUs
    # initial_M: optional, default: the activation_energy one
  microkinetics:
    initial_M: {CO: 0.05, H2: 0.05, ete: 0.05, I1_0L: 5e-4}
    simulation_time_s: 100000
    catalyst_snapshot_h: 1         # time of the catalyst-distribution snapshot
    conversion: 0.99               # yield that defines the conversion time
    plot_species: [I7_0L, I7_1L, I1_0L, I1_1L]
```

In every `initial_M`, species not listed start at 0. Do not list the studied species (it takes each value of `studied_range_M`) or the product (it starts at 0).

### 2.6 YAML details

MicroKatc reads study files with three deliberate differences from plain YAML:

- `1e-10` is a number (standard YAML 1.1 readers treat it as text).
- Only lowercase `true` and `false` are booleans, so a species called `NO` or `ON` works as a key.
- A key given twice is an error, not silently overwritten.

## 3. Check and run

```bash
python microkatc.py check my_study/study.yaml     # checks everything, runs nothing
python microkatc.py run my_study/study.yaml       # runs every requested analysis
python microkatc.py run my_study/study.yaml --fresh   # deletes my_study/results first
```

`check` reports every problem at once, each with its location in the file. `run` does the same checks first.

| Exit code | Meaning |
| --- | --- |
| 0 | Success |
| 1 | Invalid study, or `results/` made from a different version of the study |
| 2 | Failure during the run (thermochange, COPASI) |

The hydroformylation example takes about 3.5 minutes on a workstation and 7 minutes on a 4-core machine; most of it is the degree of rate control, which simulates the network twice per step and per concentration. `python main.py` is a shortcut for running it.

**Re-running.** Results are saved and reused. `results/run_info.json` records a hash of the study file and of every energy input. If you edit either, the next `run` stops and asks you to use `--fresh`, so old results are never reused for a changed study. A run that was interrupted can simply be started again.

**One study per process.** If you drive MicroKatc from your own Python script, run each study in its own process (for example `subprocess.run([sys.executable, "microkatc.py", "run", path])`): COPASI resolves relative output paths against the first folder it used in a process.

## 4. Outputs

Everything goes to `results/` next to the study file.

| Path | Contents |
| --- | --- |
| `run_info.json` | Study and input hashes, MicroKatc and thermochange commits, package versions, start time and duration |
| `G_values_of_compounds/` | Gibbs energy of every species at each temperature (kcal mol<sup>-1</sup>) |
| `G_values_of_reactions/` | Forward and reverse barrier of every step at each temperature |
| `microkinetics_simulations/` | Raw COPASI time courses, and the `df_flux_*`, `df_rate_*` and `df_drc_*` tables with every fitted E<sub>a</sub>, R<sup>2</sup> and DRC |
| `microkinetics_simulations_images/` | Arrhenius plots (R<sup>2</sup> > 0.9), E<sub>a</sub> of the product-forming steps and of product formation, DRC of every step, catalyst distribution, concentration profiles and conversion time |
| `microkinetics.json` | Catalyst amount per cycle at the snapshot time and time to the conversion threshold, per studied concentration |
| `model.cps` | The COPASI model of the last simulation |

## 5. Troubleshooting

Problems in the study file are reported before anything runs. Examples:

| Message | Fix |
| --- | --- |
| `conditions.temprature_K: unknown key; did you mean temperature_K?` | Fix the spelling. |
| `steps[4] "I4_0L + CO I5_0L": no <=> between the two sides` | Write the step as `A + B <=> C`. |
| `no energy for TS3_0L: GaussOutputFiles/TS3_0L.out not found` | Add the file, or fix the name in the step. |
| `species.energies_kcal_mol.TS1_0L: needs values at 325, 337.5, 362.5, 375 K for activation_energy` | Give one value per temperature, or drop the activation-energy analysis. |
| `analyses.microkinetics.initial_M.PMe3: the studied species is set by conditions.studied_range_M` | Remove it from `initial_M`. |
| `analyses.activation_energy.sampling_time_h: 3 h is after the end of the simulation (10000 s = 2.8 h)` | Sample earlier or simulate longer. |
| `$thermochange is not set: export thermochange=/path/to/thermochange` | Export the variable. |

Warnings do not stop the run: species in no cycle (not counted in the catalyst distribution), negative barriers, and the overall ΔG from your energies (−25.9 kcal mol<sup>-1</sup> for the example, the paper's value), printed as a check.

During a run:

| Message | Cause and fix |
| --- | --- |
| `thermochange gave a corrected G of 0.0 hartree for ...` | thermochange's own Python step failed, usually because `numpy` or `scipy` is missing in the Python it uses. Its error output follows the message. |
| `... couldn't converge to a solution. Resimulating...` | COPASI stopped early; MicroKatc retries up to 5 times. |
| `ConvergenceError: Solution can't be found for ...` | All retries failed: shorten the simulation or use less extreme concentrations. |
| `C(...) = ...M couldn't reach the threshold ...` | The product never reached the conversion threshold in the simulated time; that point is left out of the conversion-time figure. Simulate longer. |
| `Skipping Ea fit of ...: sign changes across the temperature range` | That step's net flux changes direction between temperatures, so it has no meaningful Arrhenius slope. Expected for near-equilibrium steps. |

Run the test suite with `for t in tests/test_*.py; do python "$t"; done`. `tests/test_paper_barriers.py` needs `$thermochange` and is skipped without it.
```

- [ ] **Step 2: Apply the README and CLAUDE.md patch**

Save the block below exactly as `task10.patch` in the repository root, then apply it:

```diff
--- a/README.md
+++ b/README.md
@@ -40,8 +40,9 @@
 
 | Module | Role |
 | --- | --- |
-| [`main.py`](main.py) | Entry point; all analysis parameters are set here |
-| [`get_G_compounds.sh`](get_G_compounds.sh), [`calculating_G_for_microkinetics.py`](calculating_G_for_microkinetics.py) | Thermochemistry: G of each species and the forward/reverse barrier of each step |
+| [`microkatc.py`](microkatc.py) | The command: `run` and `check` a study file |
+| [`study.py`](study.py), [`steps.py`](steps.py), [`study_yaml.py`](study_yaml.py) | Read and check the study file |
+| [`thermochemistry.py`](thermochemistry.py), [`calculating_G_for_microkinetics.py`](calculating_G_for_microkinetics.py) | G of each species (thermochange or typed) and the forward/reverse barrier of each step |
 | [`microkinetics_simulation.py`](microkinetics_simulation.py) | COPASI simulations with convergence checks and caching; catalyst and conversion analyses |
 | [`apparent_activation_energy.py`](apparent_activation_energy.py) | Apparent E<sub>a</sub> fits and DRC calculation |
 | [`plotting_functions.py`](plotting_functions.py) | All figures |
@@ -60,7 +61,7 @@
 - The rate-determining step moves with the conditions: I8_0L ⇌ I9_0L controls the rate at low PMe<sub>3</sub>, and I3_1L ⇌ I4_1L takes over at high PMe<sub>3</sub>. The negative DRC of I3_0L ⇌ I4_0L shows that this step inhibits the 0L cycle.
 - More PMe<sub>3</sub> also releases more CO, which poisons the catalyst at the I6 ⇌ I7 steps. This is why the E<sub>a</sub> of the rate-determining steps rises and then plateaus.
 
-All figures below are built by [`readme_figures.py`](readme_figures.py) from the results `main.py` saves. Blue is the 0L cycle, pink the 1L cycle and green the product. `main.py` also writes the full versions, for every step and every concentration, to `microkinetics_simulations_images/`.
+All figures below are built by [`readme_figures.py`](readme_figures.py) from the results of the example study. Blue is the 0L cycle, pink the 1L cycle and green the product. Every run also writes the full versions, for every step and every concentration, to `results/microkinetics_simulations_images/`.
 
 ### 1. Arrhenius check
 
@@ -68,7 +69,7 @@
   <img width="100%" alt="ln(r) against 1000/T for the rate-determining step and product formation, in the 0L and 1L regimes" src="pics/arrhenius.png"/>
 </p>
 
-The apparent E<sub>a</sub> comes from the slope of ln(r) against 1/T over 325–375 K, which is valid at low catalyst concentration. In each regime, product formation follows the same straight line as its rate-determining step, so the two share one E<sub>a</sub>. `main.py` draws this plot for every step and keeps those with R<sup>2</sup> > 0.9.
+The apparent E<sub>a</sub> comes from the slope of ln(r) against 1/T over 325–375 K, which is valid at low catalyst concentration. In each regime, product formation follows the same straight line as its rate-determining step, so the two share one E<sub>a</sub>. Every run draws this plot for every step and keeps those with R<sup>2</sup> > 0.9.
 
 ### 2. Apparent activation energy
 
@@ -120,7 +121,7 @@
 
 ## Reproducing the paper
 
-With the settings in `main.py`, a full run takes a few minutes and reproduces the published results:
+Running [`examples/hydroformylation/study.yaml`](examples/hydroformylation/study.yaml) takes a few minutes and reproduces the published results:
 
 | Published result | Paper | This code |
 | --- | --- | --- |
@@ -129,7 +130,7 @@
 | Apparent E<sub>a</sub> of product formation (Figure 6) | 23.5 → 21.4 → 22.1 kcal mol<sup>-1</sup> | 23.5 → 21.5 → 22.1 kcal mol<sup>-1</sup> |
 | Minimum DRC of I3_0L ⇌ I4_0L (Figure 5) | −0.22 at 3 × 10<sup>-4</sup> M | −0.22 at 3 × 10<sup>-4</sup> M |
 
-[`tests/test_paper_barriers.py`](tests/test_paper_barriers.py) checks the barriers against SI Table S3 on every run with thermochange installed.
+[`tests/test_paper_barriers.py`](tests/test_paper_barriers.py) checks the barriers against SI Table S3 on every pull request, and a nightly run checks the full results against the paper.
 
 ## Quick start
 
@@ -140,10 +141,10 @@
 python3.10 -m venv .venv && source .venv/bin/activate
 pip install -r requirements.txt                          # includes COPASI's Python bindings and copasi_helper
 export thermochange=/path/to/thermochange
-python main.py                                           # runs the hydroformylation example of the paper
+python microkatc.py run examples/hydroformylation/study.yaml   # the paper's example
 ```
 
-**To study your own reaction,** see the **[usage guide](docs/USAGE.md)**. It covers the input files and naming rules, every parameter in `main.py`, the outputs, re-running after a change, and troubleshooting.
+**To study your own reaction,** describe it in one YAML study file: the steps (`A + B <=> C via TS`), where the energies come from (Gaussian output files or typed Gibbs energies), the conditions and the analyses. Then run `python microkatc.py run my_study.yaml`. The **[usage guide](docs/USAGE.md)** covers every key, both [examples](examples/), the outputs and troubleshooting.
 
 ## Modelling assumptions
 
--- a/CLAUDE.md
+++ b/CLAUDE.md
@@ -9,10 +9,12 @@
 
 ## Pipeline
 
-1. `get_G_compounds.sh T P` runs thermochange (`$thermochange` env var) on every
-   `GaussOutputFiles/*.out` and calls `calculating_G_for_microkinetics.py`.
-2. `main.py` reads `reactions.csv`, runs COPASI simulations through `copasi_helper`
-   (`microkinetics_simulation.py`), then `apparent_activation_energy.py` and `plotting_functions.py`.
+1. `microkatc.py run <study.yaml>` reads and checks the study (`study.py`, `steps.py`,
+   `study_yaml.py`; format in `docs/USAGE.md` and `docs/superpowers/specs/2026-10-02-yaml-study-input-design.md`).
+2. `thermochemistry.py` gets G per species from thermochange (`$thermochange` env var) or typed
+   values, and writes the barrier tables in `<study folder>/results/`.
+3. The analyses run there: COPASI simulations through `copasi_helper` (`microkinetics_simulation.py`),
+   then `apparent_activation_energy.py` and `plotting_functions.py`.
 
 thermochange, copasi_helper and COPASI are external and are not installed in CI. The full
 pipeline cannot run in GitHub Actions; check changes with `python -m py_compile *.py`, `ruff check`,
@@ -23,7 +25,9 @@
 - Keep the scientific results identical. Do not change formulas, constants
   (`HARTREE_TO_KCAL_MOL`, the 4 kcal/mol barrierless barrier), default parameters, units, or the
   order of operations. A refactor that could change a number is out of scope.
-- Do not edit, rename or delete files in `GaussOutputFiles/`, `reactions.csv` or `pics/`.
+- Do not edit, rename or delete files in `examples/hydroformylation/GaussOutputFiles/` or `pics/`.
+- Run one study per Python process: COPASI resolves relative output paths against the first folder
+  it used in a process.
 - Do not rename public modules or entry points without updating every caller (`grep` all `.py`
   and `.sh` files) and the README.
 - No new runtime dependencies. Dev tools (ruff) are fine.
@@ -32,7 +36,8 @@
 
 ## Reproducing the paper
 
-With thermochange exported and the pinned requirements (Python 3.10), `python main.py`
-reproduces the paper's figures, and `tests/test_paper_barriers.py` checks the barriers against SI
-Table S3. Keep `main.py`'s conditions (0.05 M reactants, [Rh] = 1e-6 M for Ea and DRC, 5e-4 M for
-the catalyst distribution, 350 K, 19 concentrations from 1e-10 to 0.1 M) unless asked to change them.
+With thermochange exported and the pinned requirements (Python 3.10),
+`python microkatc.py run examples/hydroformylation/study.yaml` reproduces the paper's figures;
+`tests/test_paper_barriers.py` checks the barriers against SI Table S3, and
+`tests/check_reproduction.py` and `tests/check_equivalence.py` check a full run. Keep the example
+study's conditions unless asked to change them.
```

```bash
git apply --check task10.patch && git apply task10.patch && rm task10.patch
grep -n "main.py\|reactions.csv\|get_G_compounds" README.md CLAUDE.md docs/USAGE.md
```
Expected: only `docs/USAGE.md` mentions `main.py`, as the shortcut for the example.

- [ ] **Step 3: Check every link in the docs resolves**

```bash
python - <<'CHECK'
import re, pathlib
for doc in ("README.md", "docs/USAGE.md", "CLAUDE.md"):
    base = pathlib.Path(doc).parent
    for target in re.findall(r"\]\(([^)#]+)\)", pathlib.Path(doc).read_text()):
        if not target.startswith("http") and not (base / target).exists():
            print(f"{doc}: broken link {target}")
print("links checked")
CHECK
```
Expected: `links checked` and nothing else.

- [ ] **Step 4: Commit**

```bash
git add docs/USAGE.md README.md CLAUDE.md
git commit -m "docs: usage guide and README for the YAML study file"
```

---

### Task 11: Full verification against the baseline and the paper

**Files:** none changed (verification only).

- [ ] **Step 1: Run the example from a clean results folder**

```bash
rm -rf examples/hydroformylation/results
MPLBACKEND=Agg python main.py; echo "exit $?"
```
Expected: `Overall reaction ete + CO + H2 <=> prod: ΔG = -25.9 kcal/mol at 350.0 K`, `Results in .../examples/hydroformylation/results`, `exit 0`, in about 3.5 minutes (7 on 4 cores).

- [ ] **Step 2: Compare with the paper and with the baseline**

```bash
python tests/check_reproduction.py | tail -1
python tests/check_equivalence.py examples/hydroformylation/results
```
Expected: `All results match the paper.` and `Run matches the baseline.` Measured in a dry run of this plan on 2026-10-02: barriers 0 difference; Ea of the rate-determining steps 3.8 × 10<sup>-5</sup>; Ea of product formation 0.0072; DRC of the key steps 5.5 × 10<sup>-4</sup> up to 10<sup>-2</sup> M and 0.021 above; catalyst amounts 1.9 × 10<sup>-6</sup> relative; times to 99 % identical.

- [ ] **Step 3: The README figures still build**

```bash
MPLBACKEND=Agg python readme_figures.py | grep -c '^Saved'
git checkout -- pics
```
Expected: `7`. The regenerated images are discarded: the figures are refreshed separately, never as a side effect of this change.

- [ ] **Step 4: Typed example, check command and lint**

```bash
python microkatc.py check examples/typed_energies/study.yaml
python microkatc.py run examples/typed_energies/study.yaml && rm -rf examples/typed_energies/results
ruff check . && ruff format --check .
git status --short
```
Expected: `OK (2 steps, 4 species)`, a successful run, ruff clean, and an empty `git status` (results folders are git-ignored).

- [ ] **Step 5: Open the pull request**

Push the branch and open one pull request against `main`. Its body lists the spec and plan paths, the Task 11 Step 2 output, and the measured tolerances. The nightly workflow runs on it (it changes `*.py`), and must pass both checks before merging.
