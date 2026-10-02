# Using MicroKatc on your own system

This guide takes you from your own DFT calculations to the full set of analyses. The repository ships a complete worked example (the Rh-catalysed hydroformylation of the paper): run it once as is before changing anything, so you know your installation reproduces the published numbers.

- [1. Install](#1-install)
- [2. Prepare the inputs](#2-prepare-the-inputs)
- [3. Set the parameters](#3-set-the-parameters)
- [4. Run](#4-run)
- [5. Outputs](#5-outputs)
- [6. Re-running after you change something](#6-re-running-after-you-change-something)
- [7. Troubleshooting](#7-troubleshooting)

## 1. Install

You need Python 3.10, Bash, and thermochange. The COPASI Python bindings and copasi_helper are installed by `requirements.txt`.

```bash
git clone https://gitlab.com/dgarayr/thermochange.git
git clone https://github.com/0rkhann/MicroKatc.git
cd MicroKatc
python3.10 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export thermochange=/path/to/thermochange
```

`thermochange` must be **exported**, not just set: `main.py` calls `get_G_compounds.sh` in a subprocess, which only sees exported variables.

Check the installation by running the example and comparing it with the paper:

```bash
python main.py
python tests/check_reproduction.py
```

The last line should read `All results match the paper.`

## 2. Prepare the inputs

MicroKatc reads two inputs from the repository root.

### `GaussOutputFiles/`

One Gaussian output file per chemical species and per transition state, from a frequency calculation (thermochange needs the vibrational frequencies to correct the Gibbs energies). **The file name, without `.out`, is the species name** used everywhere else. The free molecules (reactants, products, ligands) need their own files too.

```
GaussOutputFiles/
├── CO.out          # free molecules
├── H2.out
├── ete.out
├── PMe3.out
├── prod.out
├── I1_0L.out       # catalyst intermediates
├── ...
└── TS1_0L.out      # transition states
```

### `reactions.csv`

One row per elementary step:

| Column | Content |
| --- | --- |
| `Rx` | The step, as `A + B = C`. Species separated by ` + `, sides by one `=`. Every step is reversible. |
| `TS` | The name of the transition state's output file, or `-` for a barrierless step. |
| `Gdir`, `Ginv` | Leave empty: MicroKatc fills in the forward and reverse barriers. |

```csv
Rx,TS,Gdir,Ginv
I1_0L = I2_0L + CO,-,,
I2_0L + ete = I3_0L,-,,
I3_0L = I4_0L,TS1_0L,,
I6_0L + H2 = I8_0L,TS3_0L,,
```

A barrierless step is treated as diffusion-controlled: its transition state is placed 4 kcal mol<sup>-1</sup> above the higher of its two sides (Besora et al., 2018).

### Naming rules

The analyses find species by name, so three conventions matter:

- **The product must be called `prod`.** The DRC is computed on the rate of `prod`, and the conversion time watches its concentration.
- **Catalyst intermediates are named `I<number>[letters]_<cycle>`**, for example `I3_0L`, `I2c_1L` or `I9t_1L`. The catalyst distribution sums every intermediate whose name ends in `_<cycle>`, for each label in `cycles`.
- **Every species in `Rx` and every name in `TS`** needs a matching `.out` file.

## 3. Set the parameters

All parameters are set, and commented, at the top of `main()` in [`main.py`](../main.py). Concentrations are in M.

**Shared**

| Parameter | Meaning | Example |
| --- | --- | --- |
| `temperature_value` | Working temperature (K) | `350.0` |
| `reactant_to_study` | The species whose initial concentration is varied | `"PMe3"` |
| `reactant_concentration_array` | Its initial concentrations; set by `left_border_concentration`, `right_border_concentration` (log10 M) and `num` | 19 points, 10<sup>-10</sup>–10<sup>-1</sup> M |
| `time_step` | Output step of the simulations (s) | `1` |

**Apparent activation energy and DRC.** These need a low catalyst concentration (a differential reactor), so that rates stay in the low-conversion regime.

| Parameter | Meaning | Example |
| --- | --- | --- |
| `c0_1` | Initial concentrations: reactants, `prod`, the catalyst's starting intermediate | `{"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 1e-6, ...}` |
| `T_values_array_Ea` | Temperatures for the Arrhenius fits (K) | 325–375 K, 5 points |
| `total_simulation_time_Ea` | Simulated time (s) | `10_000` |
| `time` | Time at which rates are read for the fits (h) | `2` |
| `e_shift` | Barrier shift for the central-difference DRC (kcal mol<sup>-1</sup>) | `0.1` |
| `total_simulation_time_drc` | Simulated time for the DRC (s); rates are read at half of it | `10_000` |
| `cores` | Processes for the DRC | `8` |

**Catalyst distribution, concentration profiles and conversion time.** These use a realistic catalyst concentration.

| Parameter | Meaning | Example |
| --- | --- | --- |
| `c0_2` | Initial concentrations | `{"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 5e-4, ...}` |
| `cycles` | Cycle labels, the suffixes of the intermediate names | `["0L", "1L"]` |
| `total_simulation_time_microkinetics` | Simulated time (s) | `100_000` |
| `main_reactants` | Reactants that can limit the yield; the smallest initial concentration sets the 100 % yield | `["CO", "H2", "ete"]` |
| `compounds_to_plot` | Species in the concentration-profile figure | `["I7_0L", "I7_1L", "I1_0L", "I1_1L"]` |
| `catalyst_concentration_time` | Time of the catalyst-distribution snapshot (h) | `1` |
| `percentage_of_convertion` | Yield that defines the conversion time | `0.99` |

`reactant_to_study` must be a key of both `c0_1` and `c0_2`: its value there is replaced by each concentration of `reactant_concentration_array`.

## 4. Run

```bash
python main.py
```

The example takes about 3.5 minutes on a workstation and about 7 minutes on a 4-core machine (GitHub's CI runner). Most of that time is the DRC: it simulates the whole network twice per step and per concentration.

To rebuild the README-style summary figures after a run:

```bash
python readme_figures.py
```

## 5. Outputs

Everything is written under the repository root.

| Folder or file | Contents |
| --- | --- |
| `G_values_of_compounds/` | Gibbs energy of every species at each T and P (kcal mol<sup>-1</sup>) |
| `G_values_of_reactions/` | `reactions.csv` with `Gdir` and `Ginv` filled in, at each T and P |
| `compounds.csv` | The list of species |
| `microkinetics_simulations/` | Raw COPASI time courses (concentrations, rates, fluxes), and the `df_flux_*`, `df_rate_*` and `df_drc_*` tables with every fitted E<sub>a</sub>, R<sup>2</sup> and DRC |
| `microkinetics_simulations_images/` | All figures from `main.py`: Arrhenius plots of every step with R<sup>2</sup> > 0.9, E<sub>a</sub> of the product-forming steps and of product formation, DRC of every step, catalyst distribution, concentration profiles and conversion time |
| `pics/*.png` | The summary figures from `readme_figures.py` |

## 6. Re-running after you change something

MicroKatc saves every simulation and table, and reuses them when it finds them.

- **Parameters in `main.py`:** the table file names include a hash of the parameters that shaped them, so a change triggers a recompute on its own.
- **Inputs (`GaussOutputFiles/` or `reactions.csv`):** these are *not* part of any file name. Old Gibbs energies, simulations and tables would be reused silently. Delete the generated results first:

```bash
rm -rf G_values_of_compounds G_values_of_reactions compounds.csv microkinetics_simulations microkinetics_simulations_images
```

## 7. Troubleshooting

| Message | Cause and fix |
| --- | --- |
| `get_G_compounds.sh did not create .../G_values_...csv. Is $thermochange exported, and does every species ...` | The Gibbs-energy step failed; the script's own output follows the message. `KeyError: 'TS3_0L'` there means that species or transition state has no `TS3_0L.out` in `GaussOutputFiles/` (check the spelling). Otherwise, `export thermochange=/path/to/thermochange`, and check that every `.out` file comes from a frequency calculation. |
| `... couldn't converge to a solution. Resimulating...` | COPASI stopped before the end of the simulation; MicroKatc retries up to 5 times. |
| `ConvergenceError: Solution can't be found for ...` | All retries failed. Try a shorter `total_simulation_time_*` or less extreme concentrations. |
| `C(...) = ...M couldn't reach the threshold ...` | The product never reached `percentage_of_convertion` of the maximum yield in the simulated time; that point is left out of the conversion-time figure. Increase `total_simulation_time_microkinetics`. |
| `Skipping Ea fit of ...: sign changes across the temperature range` | That step's net flux changes direction between temperatures, so it has no meaningful Arrhenius slope. Expected for near-equilibrium steps. |
| `ValueError: COPASI returned N '.Flux' columns (or DRC coefficients) but reactions.csv has M reactions` | The saved results and `reactions.csv` disagree; delete the generated results (section 6) and run again. |

Run the test suite with `for t in tests/test_*.py; do python "$t"; done`. `tests/test_paper_barriers.py` needs `$thermochange`, and it is skipped without it.
