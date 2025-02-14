<h1 align="center">MicroKatc</h1>

<p align="center">
  <i>A fully automated tool for analyzing a system consisting of catalytic cycles using microkinetics approach</i>
  <br/><br/>
  <img width="600" alt="MesoKinetix" src="https://github.com/0rkhann/MesoKinetix/blob/main/pics/logo.png"/>
</p>

### Table of Contents

- [Quick Start](#quick-start)
    - [Install](#install)
    - [Requirements](#requirements)
    - [Run](#run)
- [Some Notes](#some-notes)
- [Example](#example)
    - [Results](#results)
    - [Additional Results](#additional-results)
- [Acknowledgements](#acknowledgements)

---

## Quick Start

### Install

To get started, you'll need to install the following packages:

1. **Thermochange**: Package for thermochemical correction
    ```bash
    git clone https://gitlab.com/dgarayr/thermochange.git
    ```

2. **COPASI Helper**: Package to handle COPASI output (microkinetics)
    ```bash
    git clone https://gitlab.com/dgarayr/copasi_helper.git
    ```
    > **_NOTE:_**  Download the COPASI software as well

3. **MicroKatc**: Package for analysis of catalytic cycles
    ```bash
    git clone https://github.com/0rkhann/MicroKatc.git
    ```

### Requirements

To install the required libraries, use:

```bash
pip install -r requirements.txt
```

> **_NOTE:_**  It is advisable to do it in the virtual environment as some packages are compatible only with older verisons of numpy and pandas

### Run
1. Place all computed .out files into the GaussOutputFiles folder
2. Create a reactions.csv file with all reactions in the cycle (refer to the example provided)
3. Define an environment variable *thermochange* pointing at the route where **thermochange** is installed:
   ```bash
   thermochange=/home/user/programs/thermochange
   ```
4. Adjust the input parameters in main.py if necessary and execute:
```bash
python3 ./main.py
```

## Some Notes

1. **Molecular Computations**: Molecules were computed using DFT with the 6-311g(d,p) basis set and ωB97XD functional.

2. **Ideal Gas Approximation**: The pressure was adjusted based on the temperature to maintain a concentration of 1M, approximating the conditions of a liquid reaction medium, while considering ideal gas behavior:
    <p align="center">
        <span style="font-size: 1.5em;">$C = \frac{P}{RT}$</span>
    </p>
    
3. **Energy Barrier**: According to Besora et al. 2018, for "barrierless" steps, steps controlled by diffusion, the energy barrier was set to 4 kcal/mol.

4. **Simplification of Cycles**: The reactions between intermediates were considered in the most straightforward manner, without accounting for potential complexities introduced by many other possible transformations between cycles or between non-neighboring intermediates.

5. **Catalyst Concentration**: For the first part (Ea, apparent activation energy), the concentration of a catalyst should be assumed to be low to substitute $\( k_i \)$ with flux/rate in the linear form of the Arrhenius equation:

<p align="center">
        <strong>Arrhenius Equation:</strong><br>
        <span style="font-size: 1.5em;">$k = A \exp\left(\frac{-E_a}{RT}\right)$</span>
    </p>
    <p align="center">
        <strong>Linear Form of Arrhenius Equation:</strong><br>
        <span style="font-size: 1.5em;">$\ln(k_i) = \ln(A) - \frac{E_a}{RT}$</span>
    </p>
    <p align="center">
        <strong>At Low Catalyst Concentration:</strong><br>
        <span style="font-size: 1.5em;">$\ln(r_i) = \ln(A) - \frac{E_a}{RT}$</span>
        </br>
        <span style="font-size: 1.5em;">$\ln(v_i) = \ln(A) - \frac{E_a}{RT}$</span>
    </p>

## Example

<p align="center">
  <img width="1000" alt="cycle_example" src="https://github.com/0rkhann/MesoKinetix/blob/main/pics/catalytic_cycle.jpg"/>
</p>

</br>

The hydroformylation catalytic cycle with a homogeneous rhodium-based catalyst was studied by varying the initial concentration of PMe3 (i.e., <b>reactant_to_study = PMe3</b>).

### Results

#### 1. ln(ri) vs. T<sup>-1</sup> at specific concentration of studied reactant

<p align="center">
  <img width="1600" alt="ln(ri)_1_T" src="https://github.com/0rkhann/MicroKatc_private/blob/main/pics/ln(ri)_1_T_PMe3_1.000000e-05_T_range_325.0K_375.0K.svg"/>
</p>

This plot assesses the linearity of different steps and checks consistency with the low catalyst concentration assumption. Only steps with an R<sup>2</sup> value greater than 0.9 are displayed by default (can be adjusted).

#### 2. Ea vs. c0(PMe3) flux based

<p align="center">
  <img width="1600" alt="Ea_c0_flux" src="https://github.com/0rkhann/MicroKatc_private/blob/main/pics/Ea_of_rate_determining_steps.png"/>
</p>

It was interesting to observe how the activation energy of the rate-determining steps (rate determining steps were identified by the DRC vs. c0 plot) changes with increasing initial concentrations of PMe3. Initially, the activation energy increases and then plateaus as the catalyst becomes saturated with PMe3. The activation energy in the 1L cycle is lower than in the 0L cycle for both types of steps, indicating that the 1L cycle is thermodynamically more favorable and, therefore, has faster kinetics, making it the predominant pathway for product formation. The increasing trend of activation energy until the plateau can be attributed to the fact that as more PMe3 is introduced into the system, more CO is released, leading to greater catalyst poisoning during the I6 ⇌ I7 steps.

#### 3. Ea vs. c0(PMe3) compound based

<p align="center">
  <img width="700" alt="Ea_c0_compound" src="https://github.com/0rkhann/MesoKinetix/blob/main/pics/Ea_c0(PMe3)_rate_based.svg"/>
</p>

It was interesting to observe how the activation energy, related to the rate of change in compound concentrations, varies with increasing initial concentrations of PMe3. Overall, the activation energy decreases as PMe3 concentration increases, indicating greater activation and thermodynamic favorability of the 1L cycle. The slight increase in activation energy at higher concentrations may be attributed to catalyst poisoning.

#### 4. DRC vs. c0(PMe3)

<p align="center">
  <img width="1000" alt="DRC_c0" src="https://github.com/0rkhann/MicroKatc_private/blob/main/pics/DRC_of_influencing_steps.png"/>
</p>

It is important to determine the rate-determining steps to thoroughly analyze the system. The Degree of Rate Control (DRC) is a technique used for this purpose. By examining a DRC plot, one can identify these key steps. Additionally, it is interesting to observe how the DRC of different steps changes with increasing initial concentrations of PMe3.

From the plot, we can see that the step I8_0L = I9_0L is a rate-determining step in the 0L cycle. However, as the concentration of PMe3 increases, this step becomes less significant for the system because the 1L cycle starts to become more active. In contrast, for the step I3_1L = I4_1L, its DRC increases with higher concentrations of PMe3, making it a more rate-determining step. Another notable observation is that the DRC of the I3_0L = I4_0L step becomes negative with increased PMe3 concentration, indicating inhibition in the 0L cycle.

#### 5. c(catalyst) vs. c0(PMe3)

<p align="center">
  <img width="700" alt="c_catalyst_c0" src="https://github.com/0rkhann/MicroKatc_private/blob/main/pics/C(catalyst)_C(PMe3)_350.0K_2.91006e%2B03atm.svg"/>
</p>

This plot shows the concentration of the catalyst at a specified time (adjusted in case) across different cycles with increasing initial concentrations of PMe3. It is evident that as more PMe3 is introduced into the system, the 1L cycle becomes more predominant. After a certain point, the concentration of the catalyst in the 1L cycle increases significantly, surpassing that of the 0L cycle.

#### 6. Concentration evolution over time

<p align="center">
  <img width="1000" alt="c_evolution_c0" src="https://github.com/0rkhann/MicroKatc_private/blob/main/pics/C(compounds)_C(PMe3)_350.0K_2.91006e%2B03atm.svg"/>
</p>

This plot shows the evolution of concentrations for different species (adjusted in case) with varying initial concentrations of PMe3. The red line, representing the product, reaches its maximum concentration more quickly with higher amounts of PMe3 introduced into the system. This indicates that the 1L cycle is more efficient at producing the product and does it faster.

#### 7. Time of product conversion vs. c(PMe3)

<p align="center">
  <img width="700" alt="t_conversion_c" src="https://github.com/0rkhann/MicroKatc_private/blob/main/pics/C(PMe3)_total_t_product_conversion_350.0K_2.91006e%2B03atm.svg"/>
</p>

This plot shows the time required for the product to reach a specified threshold of 99% (adjusted in case) with varying initial concentrations of PMe3. It can be observed that as the amount of PMe3 introduced into the system increases, the time for the product to reach 99% of its maximum concentration decreases significantly, considering the concentration of the limiting reactant. This metric indicates that the 1L cycle is more efficient at producing the product and does so more quickly.

### Additional Results

These results are not provided in the script. However, they can easily be obtained by adding a few lines of code.

<p align="center">
  <img width="1000" alt="Ea_c0" src="https://github.com/0rkhann/MicroKatc_private/blob/main/pics/final_pieces.png"/>
</p>

> **_NOTE:_**  The first plot was simulated at a low concentration of the catalyst to meet the conditions required for the apparent activation energy assumption.

This figure was created as a final piece to complement the overall analysis. The first plot shows the change in concentration of "poisoning" intermediates over time. It demonstrates that the concentrations of I1_0L and I7_0L decrease, while those of I1_1L and I7_1L increase as the initial concentration of PMe3 rises.

From the second plot, we observe that as more PMe3 is introduced into the system, the activation energy (Ea) of the rate-determining steps in both the 0L and 1L cycles increases. However, the 1L cycle becomes predominant in product release, which is why we see a decrease in the product’s Ea on the plot of Ea versus c0(PMe3). What isn’t immediately obvious from this plot is why the Ea of the rate-determining step in the 1L cycle does not decrease even as the 1L cycle becomes more active and the catalyst is more present in the 1L cycle. To understand this, we need to look at the first and third plots.

The third plot shows that as c0(PMe3) increases, the DRC of the I8_0L = I9_0L step decreases because the I3_0L = I4_0L step inhibits the 0L cycle, facilitating the shift to the 1L cycle. This, in turn, increases the DRC of the I3_1L = I4_1L step. As it becomes easier to enter the 1L cycle with increased c0(PMe3), the concentrations of poisoning intermediates in the 0L cycle decrease. However, the concentrations of poisoning intermediates in the 1L cycle increase, making it "harder" to form the product in the 1L cycle. This explains why the Ea for the I3_1L = I4_1L step increases.

## Acknowledgements
I sincerely thank my supervisor, Dr. Diego Ruiz Garay, and Principal Investigator, Prof. Carles Bo, for their invaluable guidance and support, which have greatly enriched my experience at ICIQ. I also extend my heartfelt thanks to my fellow summer research colleagues, who became dear friends during these past two months. The time we spent together made this journey truly rewarding and memorable.
