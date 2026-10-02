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
