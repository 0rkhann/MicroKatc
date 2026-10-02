# Design: one YAML study file as MicroKatc's input

Status: approved in conversation on 2026-10-02, section by section, then refined for gaps (section 8). Items marked **(proposed)** came out of the gap review and need approval.

## 0. Prerequisites

Built on `main` after #25 (nightly reproduction check), #26 (G-step error reporting) and #27 (usage guide) are merged. The equivalence baseline in section 6 is taken from that `main`.

Platform: Linux or macOS. thermochange and its formatters are Bash scripts; Windows works through WSL only.

## 1. Goal

Any computational chemist can describe a catalytic system in **one readable file** and run every MicroKatc analysis with **one command**, without editing Python or learning hidden naming rules. The hydroformylation example of the paper is rewritten in this format and must reproduce today's results.

**Success criteria**

- A new system needs only a study file and its energies (output files or typed values).
- `python microkatc.py run examples/hydroformylation/study.yaml` reproduces today's results within COPASI's measured run-to-run noise (two fresh runs of identical code, 2026-10-02): barrier tables exactly; E<sub>a</sub> of the rate-determining steps within 10<sup>-3</sup> and of product formation within 0.01 kcal mol<sup>-1</sup>; DRC of the key steps within 10<sup>-3</sup> up to 10<sup>-2</sup> M PMe<sub>3</sub> and within 0.03 above (ill-conditioned there: under 1 % of the catalyst in the 0L cycle) and for every other step; catalyst amounts within 10<sup>-4</sup> relative; times to 99 % within 10<sup>-6</sup>. Near-equilibrium steps, whose net flux is a small difference of large rates, vary by up to 4 kcal mol<sup>-1</sup> between identical runs and are not compared.
- Every input mistake listed in section 5 is reported before any calculation starts, with its location in the file.

## 2. Decision log

| # | Decision | Chosen | Rejected and why |
| --- | --- | --- | --- |
| D1 | Who writes the input | Any computational chemist | Only gTOFfee or ioChem-BD users: too narrow; cycle-loop authors: cannot express branches and cross-links |
| D2 | Where energies come from | Gaussian/ADF output files (through thermochange) **or** typed Gibbs energies | Files only: excludes ORCA, Q-Chem and literature values; energies only: drops the current automatic route |
| D3 | What the file holds | The whole study: network, energies, conditions, analyses | Network only: users would still edit `main.py`; network now and conditions later: postpones the main usability gain |
| D4 | Syntax | YAML with chemistry-style step strings | Extending overreact's `$scheme` format: its `<=>` means a fast equilibrium, which contradicts MicroKatc's reversible steps with barriers, and it needs a hand-written parser; cycles as ordered loops: cannot express branches |
| D5 | Example | The hydroformylation CSV example is replaced by `examples/hydroformylation/study.yaml` | Keeping both: two sources of truth |
| D6 | Mixing energy sources | One source per study | Mixing: file-derived (absolute, hartree-scale) and typed (usually relative) energies would give meaningless barriers |
| D7 | Typed energies at other temperatures | A single value is used at the working temperature only; the activation-energy analysis requires one value per temperature | Treating G as temperature-independent: silently wrong E<sub>a</sub> |
| D8 | Command | `python microkatc.py run study.yaml` (no installation) | Installable `microkatc` command: needs packaging; a separate step for a PyPI release |
| D9 | Stoichiometric coefficients | Supported (`2 A` or `2*A`); `overall_reaction` is required | Not supported: coefficients change barriers, rate laws and the maximum yield |
| D10 | YAML parsing | A strict loader: `1e-10` is a number, only lowercase `true`/`false` are booleans, duplicate keys are errors (section 3.6) | Plain `yaml.safe_load`: turns `1e-10` into a string, `NO` into `False`, and silently keeps the last of two duplicate keys (measured with PyYAML 6.0.3) |
| D11 | thermochange results | Each call checked: a corrected G of 0, missing, or over 0.1 hartree from the file's own G stops the run | Trusting the exit code: thermochange exits 0 and prints `0.000000` when its Python step fails (measured) |
| D12 | Stale results **(proposed)** | `results/run_info.json` records a hash of the study; a run on a `results/` folder made from a different study stops, unless `--fresh` is given | Documenting "delete results/ after edits" only: with typed energies, editing the study is the normal workflow, and stale reuse is silent |
| D13 | Vibrational correction **(proposed)** | `species.vibrational_correction: RRHO` (default, today's behaviour) or `Grimme` (thermochange's `-g` quasi-RRHO) | Not exposed: quasi-RRHO is now common for low-frequency modes; it is a one-flag pass-through |
| D14 | Validation-only command **(proposed)** | `python microkatc.py check study.yaml` runs every section-5 check and exits 0 or 1 without calculating | Only `run`: writing a study is iterative, and a full run takes minutes |

## 3. File format

### 3.1 Example: `examples/hydroformylation/study.yaml`

```yaml
# MicroKatc study: Rh-catalysed hydroformylation of ethylene with PMe3
# Abdullayev et al., ACS Catal. 2025, 15, 4739
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
  microkinetics:
    initial_M: {CO: 0.05, H2: 0.05, ete: 0.05, I1_0L: 5e-4}
    simulation_time_s: 100000
    catalyst_snapshot_h: 1
    conversion: 0.99
    plot_species: [I7_0L, I7_1L, I1_0L, I1_1L]
```

### 3.2 Key reference

Units are part of the key name (`_K`, `_M`, `_s`, `_h`, `_kcal_mol`); values are plain numbers.

| Key | Required | Meaning |
| --- | --- | --- |
| `name` | no | Label used in logs and figure titles |
| `species.files` | one of the two | Folder, relative to the study file, holding `<name>.out` for every species and transition state |
| `species.energies_kcal_mol` | one of the two | Typed Gibbs energies (section 3.4) |
| `species.product` | yes | The species whose formation the DRC and the conversion time follow |
| `species.overall_reaction` | yes | The net reaction, with coefficients; sets the maximum yield (section 3.5) |
| `species.cycles` | for `microkinetics` | Cycle label → list of the catalyst intermediates in that cycle |
| `species.vibrational_correction` | no | `RRHO` (default) or `Grimme` **(proposed, D13)**; files only |
| `steps` | yes | Elementary steps (section 3.3) |
| `conditions.temperature_K` | yes | Working temperature |
| `conditions.studied_species` | yes | The species whose initial concentration is varied |
| `conditions.studied_range_M` | yes | `{from, to, points}`, log-spaced; `points` ≥ 1 |
| `conditions.output_step_s` | no | Output step of every simulation, default `1` (today's `time_step`). Every sampling and snapshot time must fall on an output point (section 5) |
| `analyses.activation_energy` | no | `initial_M`, `temperatures_K {from, to, points}` (linear, at least 3), `simulation_time_s`, `sampling_time_h` |
| `analyses.degree_of_rate_control` | no | `e_shift_kcal_mol`, `simulation_time_s`, optional `initial_M` (default: the one of `activation_energy`; one of the two is required), optional `cores` (default: the number of CPUs). Runs at the working temperature; the DRC target is the product's rate at half the simulation time, as today |
| `analyses.microkinetics` | no | `initial_M`, `simulation_time_s`, `catalyst_snapshot_h`, `conversion`, `plot_species` |

At least one analysis is required; analyses left out are skipped. In every `initial_M`, species not listed start at 0, and the studied species takes each value of `studied_range_M` in turn. Neither the studied species nor the product may be listed in an `initial_M`.

**Catalyst species** are the members of `species.cycles`; the catalyst distribution sums them per cycle. A species that is in no cycle, not in `overall_reaction` and not the studied species (an off-cycle intermediate or a spectator ligand) is listed in a warning, because the distribution does not count it.

### 3.3 Steps

```
A + B <=> C                  barrierless (diffusion-controlled, 4 kcal/mol above the higher side)
A + B <=> C  via TS1         transition state from TS1's energy
2 A <=> A2  via TS2          coefficients: "2 A" or "2*A", positive integers
```

- `<=>`, `=` and `⇌` are accepted and mean the same: every step is reversible.
- Species are separated by ` + `. Names are case-sensitive, start with a letter and contain no spaces.
- A coefficient n multiplies the species' Gibbs energy in the barrier (G(TS) − n·G(A)), and the species is passed to COPASI n times. A probe on 2 A ⇌ B matched the analytical second-order solution to 8 × 10<sup>-7</sup> relative, with exact mass balance, for both `A + A` and `2*A` (2026-10-02).
- Steps are numbered in file order; that order is the one COPASI uses (`rNN`).

### 3.4 Typed energies

```yaml
species:
  energies_kcal_mol:
    I1_0L: 0.0                       # used at the working temperature only
    TS1_0L: {325: 11.5, 337.5: 11.5, 350: 11.5, 362.5: 11.5, 375: 11.5}
```

- Gibbs energies at 1 M, on any common reference (only differences enter the barriers).
- A single number is valid only at the working temperature. If `activation_energy` is requested, every species needs a value at every temperature of `temperatures_K` (matched to 10<sup>-6</sup> K); otherwise the run stops with the missing (species, temperature) pairs.
- `files` and `energies_kcal_mol` cannot be combined (D6).

### 3.5 Overall reaction and maximum yield

`overall_reaction` uses the step syntax without `via`. For a study's `initial_M`, the maximum product concentration is

min over overall reactants r of (c<sub>0</sub>(r) / ν<sub>r</sub>) × ν<sub>product</sub>,

and the conversion time is the first time the product reaches `conversion` × that value. For the example (all ν = 1) this equals today's "smallest initial reactant concentration". MicroKatc prints the overall ΔG from the energies at the working temperature as a check (−25.9 kcal mol<sup>-1</sup> in the paper).

### 3.6 YAML parsing rules (D10)

The loader is `yaml.SafeLoader` with three changes, each tested:

- **Numbers:** `1e-10` and `5e-4` are numbers. YAML 1.1, which PyYAML follows, needs a decimal point and would read them as text.
- **Booleans:** only lowercase `true` and `false`. `yes`, `no`, `on`, `off`, `Y`, `N`, `NO`, `ON` stay text, so a species called `NO` works as a key.
- **Duplicate keys** are an error naming the key and both lines. PyYAML keeps the last one silently.

Temperature keys in `energies_kcal_mol` (`325:`, `337.5:`) are read as numbers.

## 4. Architecture

```
study.yaml ─► study.py: parse + validate ─► Study object
                                               │
        ┌──────────────────────────────────────┼─────────────────────────────┐
        ▼                                      ▼                             ▼
thermochemistry.py: G per species       results/ next to the       existing analyses
(thermochange per file, or typed)       study file                 (Ea, DRC, microkinetics)
  └─► barrier table, same CSV format copasi_helper reads today
```

**New**

- `study.py`: `load_study(path) -> Study` (dataclasses). Parses YAML (PyYAML, added to `requirements.txt`) and step strings, validates (section 5), derives defaults. The only module that knows the file format.
- `thermochemistry.py`: G for every species named in the study, at each needed temperature, from thermochange or from typed energies; then the barrier table through the existing `GibbsEnergyCalculator`, with coefficients. Replaces `get_G_compounds.sh` and the command-line part of `calculating_G_for_microkinetics.py`. For files it reproduces today's call exactly:
  - pressure P = (1 M)·R·T in atm at each temperature, so G is on the 1 M standard state;
  - `formatted_energy_outputter.sh <file> <T> <P>` (plus `-g` for Grimme, D13); the 4th tab-separated field is the corrected G in hartree, × 627.509 for kcal mol<sup>-1</sup>;
  - each call runs in its own temporary folder, because the script writes `temp_summary.temp` to the current folder and does not always delete it; calls may then run in parallel;
  - checks per call (D11): the corrected G must be present, non-zero and within 0.1 hartree of the file's own G (the 3rd field); otherwise the run stops, naming the species and quoting thermochange's error output;
  - before any call: `$thermochange` must be set and contain `formatters/formatted_energy_outputter.sh`.
- `microkatc.py`:
  - `python microkatc.py run <study.yaml> [--fresh]`: loads the study, writes results to `<study folder>/results/`, runs the requested analyses and saves all figures. `--fresh` deletes `results/` first.
  - `python microkatc.py check <study.yaml>` **(proposed, D14)**: validation only.
  - Exit codes: 0 success, 1 invalid study or stale results, 2 failure during the run.
  - One study per process: COPASI resolves relative output paths against the first folder it saved a model in during a process (measured: a second study in the same process wrote into the first study's `results/`). Tests and scripts run each study in its own process.
  - `load_study(path, check_energy_sources=False)` skips the output-file and `$thermochange` checks, for scripts that only read a finished run (`readme_figures.py`, `tests/check_reproduction.py`).
  - `results/` also holds `model.cps`, the COPASI model of the last simulation, which the old code wrote to the repository root.
  - Writes `results/run_info.json`: a SHA-256 of the study file and of every energy input (output files or typed values), the MicroKatc git commit, the thermochange commit when it is a git checkout, the versions of Python, python-copasi, numpy and scipy, and the start time and duration. With D12, a `run` on a `results/` whose recorded study hash differs stops with: `results/ was made from a different version of study.yaml; delete it or run with --fresh`.

**Changed**

- `CURRENT_DIRECTORY` (fixed at import) is replaced by an explicit run directory passed to the analyses; scripts are located from the code's own folder. MicroKatc then runs from any folder.
- The product name (today `"prod"` in the DRC target, the conversion time and the plot selection) comes from `species.product`. The generic figures keep today's selection: Arrhenius plots of every step at the lowest studied concentration, E<sub>a</sub> of the steps that form the product and of product formation, DRC of every step.
- Cycle membership comes from `species.cycles` instead of the name regex.
- The conversion threshold uses section 3.5.
- Figure grids (today the fixed `nrows`/`ncols` in `main.py`) are sized from the number of panels.
- `main.py` becomes a one-line wrapper that runs the example study.
- `readme_figures.py` and `tests/check_reproduction.py` stay specific to the hydroformylation example (which steps to highlight, the 0L/1L colours, the paper's values); they read concentrations, times and paths from the example study and its `results/` instead of their own copies.

**Moved or removed**

- `GaussOutputFiles/` → `examples/hydroformylation/GaussOutputFiles/`; `reactions.csv` removed (its steps are in the study file).
- `get_G_compounds.sh` and `bash_parsing.py` removed once `thermochemistry.py` replaces them.
- `results/` folders are git-ignored.

## 5. Validation and errors

All problems are collected and reported together before any calculation, each with its location.

| Check | Example message |
| --- | --- |
| YAML syntax | `study.yaml line 14: mapping values are not allowed here` |
| Unknown or misspelled keys, with the closest valid key | `conditions.temprature_K: unknown key; did you mean temperature_K?` |
| Required keys; `files` xor `energies_kcal_mol`; at least one analysis | `species: give files or energies_kcal_mol, not both` |
| Step syntax | `steps[4] "I4_0L + CO I5_0L": no <=> between the two sides` |
| Coefficient syntax | `steps[7] "1.5 A <=> B": coefficients must be positive integers` |
| An energy source for every species and TS | `no energy for TS3_0L: GaussOutputFiles/TS3_0L.out not found` |
| Name consistency: product, studied species, cycle members, `initial_M` keys, `overall_reaction` species and `plot_species` appear in steps; cycles do not overlap; each `initial_M` contains a catalyst intermediate when `cycles` is given | `cycles.1L: I9x_1L does not appear in any step` |
| Typed energies cover the Ea temperatures | `energies_kcal_mol.TS1_0L: needs values at 325, 337.5, 362.5, 375 K for activation_energy` |
| Numbers | `analyses.activation_energy.sampling_time_h: 3 h is after the end of the simulation (10000 s = 2.8 h)`; also positive times and concentrations, `from < to`, at least 3 Ea temperatures, `0 < conversion < 1` |
| `$thermochange` unset or wrong, when `files` is used | `$thermochange is not set: export thermochange=/path/to/thermochange` |
| Studied species or product listed in an `initial_M` | `analyses.microkinetics.initial_M.PMe3: the studied species is set by conditions.studied_range_M` |
| Duplicate steps; duplicate YAML keys | `steps[9] repeats steps[2] (I9_0L <=> prod + I2_0L)`; `species.energies_kcal_mol: key TS1_0L appears twice (lines 12 and 31)` |
| Sampling and snapshot times off the output grid (the code reads rows by exact time) | `analyses.microkinetics.catalyst_snapshot_h: 1 h is not a multiple of output_step_s (7 s)` |
| Warnings (the run continues): catalyst species in no cycle; a transition state below one of its sides (negative barrier); the overall ΔG | `warning: steps[3] TS1_0L is 2.1 kcal/mol below I4_0L (negative reverse barrier)` |

Run-time errors (thermochange failure, COPASI convergence) keep their current messages and are prefixed with the study path.

## 6. Testing

1. **Equivalence (before any code change):** run `main.py` on current `main`, save the barrier tables at the 5 temperatures and the flux, rate and DRC tables. After the change, the example study must match within the section 1 tolerances, checked by `tests/check_equivalence.py <results>`. `tests/test_paper_barriers.py` is switched to build its barrier table from the study file and still checks it against SI Table S3; the nightly checks the full run.
2. **Parser and validation:** unit tests for valid step strings (`<=>`, `=`, `⇌`, `via`, coefficients, spacing) and one test per error row in section 5, asserting the message.
3. **Typed energies:** `examples/typed_energies/`, a three-step cycle with typed values; its barrier table is checked against hand-computed values.
4. **Stoichiometry:** `tests/test_stoichiometry.py` (the 2 A ⇌ B probe against the analytical solution; runs in CI, which installs COPASI), a barrier test for G(TS) − 2·G(A), and the maximum yield of 2 A + B → P.
5. **Runs from anywhere:** start the command from another folder and check inputs are found and results land in the study's `results/`.
6. **YAML loader (D10):** `1e-10` and `5e-4` load as numbers, `NO` as text, a duplicate key raises.
7. **thermochange checks (D11):** a fake thermochange that prints `0.000000` stops the run with the species named; `$thermochange` unset is reported before any call.
8. **Provenance and stale results (D12):** `run_info.json` is written with every field; a second run after editing the study stops; `--fresh` clears and runs.
9. **Existing tests** updated to the new constructors and paths: `test_cache_names.py`, `test_convergence_retry.py`, `test_drc_central_difference.py`, `test_paper_barriers.py`, and `test_g_step_failure.py` from #26 (its fake `get_G_compounds.sh` becomes a fake thermochange). Nightly runs the YAML example, and `check_reproduction.py` must still pass all 11 checks.
10. **Docs:** `docs/USAGE.md` (from #27) rewritten around the study file: key reference, YAML rules, both examples, error and warning table. The README Quick start becomes the single `microkatc.py run` command.

## 7. Out of scope

- Reading ORCA, Q-Chem or VASP output files (thermochange supports Gaussian and ADF).
- A PyPI package and an installable `microkatc` command (D8).
- Finer-grained cache invalidation: with D12, any change to the study or its energy inputs stops reuse of the whole `results/` folder; reusing the unaffected parts is not attempted.
- Irreversible steps (`->`).

## 8. Gap review (2026-10-02)

Each gap was checked against the code or measured; the resolutions are in the sections above.

| Gap | Evidence | Resolution |
| --- | --- | --- |
| `1e-10` and `5e-4` in the example parse as text | PyYAML 6.0.3: `{from: 1e-10}` → `'1e-10'` | Strict loader (3.6, D10) |
| A species named `NO` or `ON` becomes a boolean | `{NO: 0.05}` → `{False: 0.05}` | Strict loader (3.6, D10) |
| Duplicate keys are silently overwritten | `{A: 1, A: 2}` → `{'A': 2}` | Error (3.6, section 5) |
| A broken thermochange gives G = 0 for every species, with exit code 0 | numpy missing in thermochange's Python: corrected G printed as `0.000000` | Per-call checks (4, D11) |
| thermochange leaves `temp_summary.temp` in the current folder | the cleanup line is commented out in its Gaussian branch | One temporary folder per call (4) |
| The vibrational correction is not selectable | thermochange's `-g` flag | `vibrational_correction` (D13, proposed) |
| The output time step had no key | `main.py` sets `time_step = 1` for every simulation | `conditions.output_step_s` (3.2) |
| Sampling times off the output grid crash | rows are read with `df["time"] == time` | Validation (section 5) |
| The DRC could not run without the activation-energy analysis | it borrowed that analysis' `initial_M` | Its own optional `initial_M` (3.2) |
| "Catalyst species" was undefined | the distribution needs it; a spectator ligand must not be counted as catalyst | Defined as the `cycles` members; other unassigned species warned about (3.2) |
| Edited studies reuse stale results silently | caches are keyed by conditions only | `run_info.json` hash check (4, D12, proposed) |
| No record of how a result was made | — | `run_info.json` (4) |
| No way to validate without a full run | — | `check` command (D14, proposed) |
| `readme_figures.py` and `check_reproduction.py` hard-code the example's steps | `RDS_0L`, `INHIBITOR`, colours | Stay example-specific; read paths and conditions from the study (4) |
| Five existing tests depend on the constructors and paths that change | `grep` of `tests/` | Listed for update (6.9) |
| The spec depends on three open PRs | #25, #26, #27 | Prerequisites (0) |
| Negative barriers pass silently | no check in `GibbsEnergyCalculator` | Warning (section 5) |
| Equivalence tolerances were tighter than COPASI's own run-to-run noise | product E<sub>a</sub> moved 0.0044–0.0072 kcal mol<sup>-1</sup> and one DRC point 0.021 between runs whose inputs are numerically equivalent | Measured tolerances (section 1), `tests/check_equivalence.py` |
| A second study in one Python process writes into the first study's folder | reproduced while testing `microkatc.py` | One study per process (4) |
