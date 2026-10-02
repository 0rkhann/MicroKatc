# MicroKatc

Microkinetic analysis of catalytic cycles (apparent activation energy, degree of rate control,
concentration evolution) built on DFT Gibbs energies from Gaussian outputs. The code was written
for an ICIQ summer research project and later used in:

> O. Abdullayev, D. Garay-Ruiz, B. Bori-Bru, C. Bo. "Microkinetic Assessment of Ligand-Exchanging
> Catalytic Cycles." *ACS Catal.* **2025**, 15 (6), 4739–4745. https://doi.org/10.1021/acscatal.5c00348

## Pipeline

1. `microkatc.py run <study.yaml>` reads and checks the study (`study.py`, `steps.py`,
   `study_yaml.py`; format in `docs/USAGE.md`).
2. `thermochemistry.py` gets G per species from thermochange (`$thermochange` env var) or typed
   values, and writes the barrier tables in `<study folder>/results/`.
3. The analyses run there: COPASI simulations through `copasi_helper` (`microkinetics_simulation.py`),
   then `apparent_activation_energy.py` and `plotting_functions.py`.

thermochange, copasi_helper and COPASI are external and are not installed in CI. The full
pipeline cannot run in GitHub Actions; check changes with `python -m py_compile *.py`, `ruff check`,
and by reading the code paths.

## Rules for changes

- Keep the scientific results identical. Do not change formulas, constants
  (`HARTREE_TO_KCAL_MOL`, the 4 kcal/mol barrierless barrier), default parameters, units, or the
  order of operations. A refactor that could change a number is out of scope.
- Do not edit, rename or delete files in `examples/hydroformylation/GaussOutputFiles/` or `pics/`.
- Run one study per Python process: COPASI resolves relative output paths against the first folder
  it used in a process.
- Do not rename public modules or entry points without updating every caller (`grep` all `.py`
  and `.sh` files) and the README.
- No new runtime dependencies. Dev tools (ruff) are fine.
- Keep the README's scientific discussion; fix its structure, links and wording only.
- Small, reviewable PRs: one theme per PR.

## Reproducing the paper

With thermochange exported and the pinned requirements (Python 3.10),
`python microkatc.py run examples/hydroformylation/study.yaml` reproduces the paper's figures;
`tests/test_paper_barriers.py` checks the barriers against SI Table S3, and
`tests/check_reproduction.py` and `tests/check_equivalence.py` check a full run. Keep the example
study's conditions unless asked to change them.
