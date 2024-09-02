<h1 align="center">MicroKatc</h1>

<p align="center">
  <i>A fully automated tool for analyzing a system consisting of catalytic cycles using microkinetics approach</i>
  <br/><br/>
  <img width="600" alt="MesoKinetix" src="https://github.com/0rkhann/MicroKatc/blob/main/pics/logo.png"/>
</p>

### Table of Contents

- [Quick Start](#quick-start)
    - [Install](#install)
    - [Requirements](#requirements)
    - [Run](#run)
- [Some Notes](#some-notes)
- [Example](#example)
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
Due to privacy reasons, the real example is currently **unavailable** as these results will be used in a research paper. It will be uploaded once the research paper is published.

<p align="center">
  <img width="1000" alt="cycle_example" src="https://github.com/0rkhann/MicroKatc/blob/main/pics/toy_cycle.png"/>
</p>

If you run the code with the provided toy cycle, you won't get any results because the GaussOutputFiles contain randomly generated molecules, causing the microkinetics solution to fail to converge due to the lack of chemical relevance. I didn't have time to compute molecules for a different catalytic cycle. However, this example still demonstrates the workflow, which you can use in your study.

## Acknowledgements
I sincerely thank my supervisor, Dr. Diego Ruiz Garay, and Principal Investigator, Prof. Carles Bo, for their invaluable guidance and support, which have greatly enriched my experience at ICIQ. I also extend my heartfelt thanks to my fellow summer research colleagues, who became dear friends during these past two months. The time we spent together made this journey truly rewarding and memorable.
