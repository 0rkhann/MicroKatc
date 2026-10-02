<h1 align="center">MicroKatc</h1>

<p align="center">
  <b>Automated microkinetic analysis of interconnected catalytic cycles, from DFT outputs to mechanistic insight.</b>
  <br/><br/>
  <a href="https://doi.org/10.1021/acscatal.5c00348"><img alt="Published in ACS Catalysis" src="https://img.shields.io/badge/ACS%20Catal.-2025%2C%2015%2C%204739-1f6feb"/></a>
  <a href="https://doi.org/10.1021/acscatal.5c00348"><img alt="DOI" src="https://img.shields.io/badge/DOI-10.1021%2Facscatal.5c00348-blue"/></a>
  <img alt="Python 3" src="https://img.shields.io/badge/python-3-3776AB?logo=python&logoColor=white"/>
  <a href="LICENSE"><img alt="MIT License" src="https://img.shields.io/badge/license-MIT-green"/></a>
  <br/><br/>
  <img width="600" alt="MicroKatc logo" src="pics/logo.png"/>
</p>

MicroKatc turns a set of Gaussian calculations and a list of elementary steps into a complete kinetic picture of a catalytic system. It computes thermally corrected Gibbs energies, builds and runs COPASI microkinetic models, and then answers the questions a mechanistic study needs: **which steps control the rate, how the apparent activation energy responds to reaction conditions, and where the catalyst actually spends its time.**

The approach is published in:

> O. Abdullayev, D. Garay-Ruiz, B. Bori-Bru, C. Bo. "Microkinetic Assessment of Ligand-Exchanging Catalytic Cycles." *ACS Catal.* **2025**, *15* (6), 4739–4745. [doi:10.1021/acscatal.5c00348](https://doi.org/10.1021/acscatal.5c00348)

## Highlights

- **End-to-end pipeline:** raw DFT output files in, publication-ready figures out, driven by one script.
- **Apparent activation energy (E<sub>a</sub>)** from Arrhenius fits of every step's flux and every species' rate of change, swept over temperature and over the concentration of any chosen reactant.
- **Degree of rate control (DRC)** for every elementary step, by central finite differences parallelised across CPU cores, to identify the rate-determining and inhibiting steps.
- **Multi-cycle catalyst tracking:** the catalyst concentration in each competing cycle (for example 0-ligand vs. 1-ligand cycles) as conditions change.
- **Conversion-time analysis:** time to reach a set product yield as a function of reactant concentration.
- **Result caching:** simulations and fitted parameters are saved as CSV files and reused, so long parameter sweeps never repeat finished work.

## How it works

```mermaid
flowchart LR
    A["Gaussian .out files"] --> B["thermochange<br/>G(T, P) of each species"]
    C["reactions.csv<br/>elementary steps + TS"] --> D
    B --> D["Forward / reverse<br/>Gibbs barriers"]
    D --> E["COPASI model<br/>(copasi_helper)"]
    E --> F["Time-course simulations<br/>over T and c0"]
    F --> G["Apparent Ea"]
    F --> H["Degree of rate control"]
    F --> I["Catalyst distribution<br/>and conversion time"]
```

| Module | Role |
| --- | --- |
| [`main.py`](main.py) | Entry point; all analysis parameters are set here |
| [`get_G_compounds.sh`](get_G_compounds.sh), [`calculating_G_for_microkinetics.py`](calculating_G_for_microkinetics.py) | Thermochemistry: G of each species and the forward/reverse barrier of each step |
| [`microkinetics_simulation.py`](microkinetics_simulation.py) | COPASI simulations with convergence checks and caching; catalyst and conversion analyses |
| [`apparent_activation_energy.py`](apparent_activation_energy.py) | Apparent E<sub>a</sub> fits and DRC calculation |
| [`plotting_functions.py`](plotting_functions.py) | All figures |

## Case study: Rh-catalysed hydroformylation

<p align="center">
  <img width="1000" alt="Hydroformylation catalytic cycles" src="pics/catalytic_cycle.jpg"/>
</p>

The example data in this repository model the hydroformylation of ethene by a homogeneous rhodium catalyst. Two cycles compete: one without (0L) and one with (1L) a coordinated PMe<sub>3</sub> ligand. The analysis varies the initial PMe<sub>3</sub> concentration (`reactant_to_study = "PMe3"`).

**Key findings**

- Adding PMe<sub>3</sub> shifts the catalyst from the 0L cycle into the 1L cycle. The 1L cycle has lower activation energies, becomes the main source of product, and reaches 99 % conversion much faster.
- The rate-determining step moves with the conditions: I8_0L ⇌ I9_0L controls the rate at low PMe<sub>3</sub>, and I3_1L ⇌ I4_1L takes over at high PMe<sub>3</sub>. The negative DRC of I3_0L ⇌ I4_0L shows that this step inhibits the 0L cycle.
- More PMe<sub>3</sub> also releases more CO, which poisons the catalyst at the I6 ⇌ I7 steps. This is why the E<sub>a</sub> of the rate-determining steps rises and then plateaus.

### 1. ln(r<sub>i</sub>) vs. 1/T

<p align="center">
  <img width="1600" alt="ln(ri) vs 1/T" src="pics/ln_ri_vs_1_over_T.svg"/>
</p>

Arrhenius plots for each step at a fixed PMe<sub>3</sub> concentration. They check that every step is linear in 1/T and that the low-catalyst-concentration assumption holds. By default, only steps with R<sup>2</sup> > 0.9 are shown.

### 2. Apparent E<sub>a</sub> of the rate-determining steps (flux-based)

<p align="center">
  <img width="1600" alt="Ea vs c0 flux based" src="pics/Ea_of_rate_determining_steps.png"/>
</p>

The rate-determining steps are those identified by the DRC analysis (figure 4). Their activation energy rises with c<sub>0</sub>(PMe<sub>3</sub>) and then plateaus as the catalyst saturates with PMe<sub>3</sub>. For both kinds of step, E<sub>a</sub> is lower in the 1L cycle than in the 0L cycle, so the 1L cycle has faster kinetics and dominates product formation. The rise before the plateau comes from CO release: more PMe<sub>3</sub> frees more CO, which poisons the catalyst at the I6 ⇌ I7 steps.

### 3. Apparent E<sub>a</sub> of product formation (rate-based)

<p align="center">
  <img width="700" alt="Ea vs c0 rate based" src="pics/Ea_c0(PMe3)_rate_based.svg"/>
</p>

This E<sub>a</sub> comes from the rate of change of the product concentration. Overall it falls as PMe<sub>3</sub> increases, reflecting the activation of the 1L cycle. The slight increase at the highest concentrations is attributed to catalyst poisoning.

### 4. Degree of rate control vs. c<sub>0</sub>(PMe<sub>3</sub>)

<p align="center">
  <img width="1000" alt="DRC vs c0" src="pics/DRC_of_influencing_steps.png"/>
</p>

I8_0L ⇌ I9_0L is rate-determining in the 0L cycle, but it loses importance as the 1L cycle becomes active. The DRC of I3_1L ⇌ I4_1L grows with PMe<sub>3</sub>, making it the controlling step at high concentration. The DRC of I3_0L ⇌ I4_0L becomes negative, which shows that this step inhibits the 0L cycle.

### 5. Catalyst distribution between cycles

<p align="center">
  <img width="700" alt="Catalyst concentration vs c0" src="pics/catalyst_distribution_350K.svg"/>
</p>

The catalyst concentration in each cycle at t = 1 h, as PMe<sub>3</sub> increases. Beyond about 5 × 10<sup>-4</sup> M, the 1L cycle holds more catalyst than the 0L cycle.

### 6. Concentration profiles over time

<p align="center">
  <img width="1000" alt="Concentration evolution" src="pics/concentration_evolution_350K.svg"/>
</p>

The time evolution of the resting states I1 and I7 of both cycles at each PMe<sub>3</sub> concentration. As PMe<sub>3</sub> increases, the catalyst moves from the 0L species into their 1L counterparts.

### 7. Time to 99 % conversion

<p align="center">
  <img width="700" alt="Conversion time vs c0" src="pics/conversion_time_350K.svg"/>
</p>

The time for the product to reach 99 % of the maximum allowed by the limiting reactant (the threshold is adjustable). It drops from 5.7 h with the 0L cycle alone to 1.7 h once the 1L cycle takes over.

### Putting it together

<p align="center">
  <img width="1000" alt="Combined analysis" src="pics/final_pieces.png"/>
</p>

These panels are not produced by `main.py`, but a few extra lines of code generate them from the same results. The first panel was simulated at low catalyst concentration, as the apparent E<sub>a</sub> analysis requires.

- **Panel 1:** as c<sub>0</sub>(PMe<sub>3</sub>) increases, the poisoning intermediates I1_0L and I7_0L decrease, while I1_1L and I7_1L increase.
- **Panel 2:** the E<sub>a</sub> of the rate-determining steps rises in both cycles. The product's E<sub>a</sub> still falls, because the 1L cycle takes over product release.
- **Panel 3:** the DRC of I8_0L ⇌ I9_0L falls, because I3_0L ⇌ I4_0L inhibits the 0L cycle and pushes the catalyst into the 1L cycle, raising the DRC of I3_1L ⇌ I4_1L.

This explains why the E<sub>a</sub> of the 1L rate-determining step increases even though the 1L cycle becomes more active. As entering the 1L cycle gets easier, its own poisoning intermediates build up, which makes product formation in that cycle harder.

## Reproducing the paper

With the settings in `main.py`, a full run takes a few minutes and reproduces the published results:

| Published result | Paper | This code |
| --- | --- | --- |
| Gibbs barriers of all 22 steps at 5 temperatures (SI Tables S1–S5) | 0.1 kcal mol<sup>-1</sup> precision | identical |
| Time to 99 % conversion, 0L only → 1L (Figure 4) | 5.66 h → 1.69 h | 5.67 h → 1.69 h |
| Apparent E<sub>a</sub> of product formation (Figure 6) | 23.5 → 21.4 → 22.1 kcal mol<sup>-1</sup> | 23.5 → 21.5 → 22.1 kcal mol<sup>-1</sup> |
| Minimum DRC of I3_0L ⇌ I4_0L (Figure 5) | −0.22 at 3 × 10<sup>-4</sup> M | −0.22 at 3 × 10<sup>-4</sup> M |

[`tests/test_paper_barriers.py`](tests/test_paper_barriers.py) checks the barriers against SI Table S3 on every run with thermochange installed.

## Quick start

### 1. Install dependencies

MicroKatc relies on two packages by D. Garay-Ruiz:

```bash
git clone https://gitlab.com/dgarayr/thermochange.git    # thermochemical corrections
git clone https://gitlab.com/dgarayr/copasi_helper.git   # COPASI model building and simulation
git clone https://github.com/0rkhann/MicroKatc.git
```

Then, inside a Python 3.10 virtual environment, install the Python requirements. They include the COPASI Python bindings (`python-copasi`) and the pinned copasi_helper commit:

```bash
cd MicroKatc
python3.10 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Prepare inputs

1. Put the Gaussian `.out` file of every species and transition state in `GaussOutputFiles/`.
2. List the elementary steps in `reactions.csv`, one per row, with the transition state in the `TS` column. Use `-` for a barrierless step. See the provided example.

### 3. Run

```bash
export thermochange=/path/to/thermochange   # must be exported so get_G_compounds.sh can see it
python3 main.py
```

Every analysis parameter (temperatures, concentration ranges, simulation times, the reactant to study, the number of cores for DRC) is set and commented in [`main.py`](main.py).

## Modelling assumptions

1. **Electronic structure:** species were computed with DFT at the ωB97X-D/6-311G(d,p) level.
2. **Standard state:** the pressure is set from the temperature so that the ideal-gas concentration is 1 M, approximating a liquid medium: $C = \frac{P}{RT}$.
3. **Barrierless steps:** steps without a located transition state are treated as diffusion-controlled, with a barrier of 4 kcal/mol (Besora et al., 2018).
4. **Network topology:** only the listed steps between neighbouring intermediates are included. Cross-cycle or non-neighbour transformations are not modelled.
5. **Apparent E<sub>a</sub>:** at low catalyst concentration, the rate constant $k_i$ in the linearised Arrhenius equation can be replaced by the step flux $r_i$ or the species rate $v_i$:

$$\ln r_i = \ln A - \frac{E_a}{RT}, \qquad \ln v_i = \ln A - \frac{E_a}{RT}$$

## Citation

If you use MicroKatc in your work, please cite:

```bibtex
@article{Abdullayev2025,
  author  = {Abdullayev, Orkhan and Garay-Ruiz, Diego and Bori-Bru, Berta and Bo, Carles},
  title   = {Microkinetic Assessment of Ligand-Exchanging Catalytic Cycles},
  journal = {ACS Catalysis},
  year    = {2025},
  volume  = {15},
  number  = {6},
  pages   = {4739--4745},
  doi     = {10.1021/acscatal.5c00348}
}
```

Citation metadata is also in [`CITATION.cff`](CITATION.cff); GitHub shows it under "Cite this repository".

## Acknowledgements

MicroKatc was developed during a summer research stay at [ICIQ](https://iciq.org/). I sincerely thank my supervisor, Dr. Diego Garay-Ruiz, and Principal Investigator, Prof. Carles Bo, for their guidance and support. I also thank my fellow summer research colleagues, who became dear friends and made the experience truly memorable.

## License

[MIT](LICENSE)
