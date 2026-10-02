# MicroKatc

Microkinetic analysis of catalytic cycles (apparent activation energy, degree of rate control,
concentration evolution) built on DFT Gibbs energies from Gaussian outputs. The code was written
for an ICIQ summer research project and later used in:

> O. Abdullayev, D. Garay-Ruiz, B. Bori-Bru, C. Bo. "Microkinetic Assessment of Ligand-Exchanging
> Catalytic Cycles." *ACS Catal.* **2025**, 15 (6), 4739–4745. https://doi.org/10.1021/acscatal.5c00348

## Pipeline

1. `get_G_compounds.sh T P` runs thermochange (`$thermochange` env var) on every
   `GaussOutputFiles/*.out` and calls `calculating_G_for_microkinetics.py`.
2. `main.py` reads `reactions.csv`, runs COPASI simulations through `copasi_helper`
   (`microkinetics_simulation.py`), then `apparent_activation_energy.py` and `plotting_functions.py`.

thermochange, copasi_helper and COPASI are external and are not installed in CI. The full
pipeline cannot run in GitHub Actions; check changes with `python -m py_compile *.py`, `ruff check`,
and by reading the code paths.

## Rules for changes

- Keep the scientific results identical. Do not change formulas, constants
  (`HARTREE_TO_KCAL_MOL`, the 4 kcal/mol barrierless barrier), default parameters, units, or the
  order of operations. A refactor that could change a number is out of scope.
- Do not edit, rename or delete files in `GaussOutputFiles/`, `reactions.csv` or `pics/`.
- Do not rename public modules or entry points without updating every caller (`grep` all `.py`
  and `.sh` files) and the README.
- No new runtime dependencies. Dev tools (ruff) are fine.
- Keep the README's scientific discussion; fix its structure, links and wording only.
- Small, reviewable PRs: one theme per PR.

## Known issues worth fixing

- README image links point to other repos (`MesoKinetix`, `MicroKatc_private`); use relative
  `pics/...` paths.
- `requirements.txt` is a whole-system `pip freeze` (apturl, dell-recovery, ...). Reduce it to what
  the code imports: numpy, pandas, scipy, matplotlib, seaborn, and the pinned
  `copasi_helper` git dependency.
- No citation for the paper: add `CITATION.cff` and a "Citation" README section.
- Typos: `auxilary_functions.py` (rename only with every import updated), "verisons", "foe",
  "concentraiton".
