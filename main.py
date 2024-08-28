from auxilary_functions import AuxiliaryFunctions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from apparent_activation_energy import ApparentEaAnalysis, DRCAnalysis
from microkinetics_simulation import MicroKinetics

def main():
    # Load reaction data
    reaction_df = pd.read_csv("reactions.csv", sep = ",")
    reactions = AuxiliaryFunctions.reactions_number(reaction_df)

    # Specify temperature in K
    temperature_value = 350.0

    # Specify range of T for which simulations will be calculated (Ea analysis)
    T_values_array_Ea = np.linspace(start = temperature_value - 25,
                                 stop = temperature_value + 25,
                                 num = 5,
                                 endpoint = True
    )

    # Creating an array of concentrations of a reactant (lg scale)
    left_border_concentration = -10
    right_border_concentration = 0
    reactant_concentration_array = np.logspace(start = left_border_concentration, 
                                               stop = right_border_concentration, 
                                               num = 9, 
                                               endpoint = True
    )

    # Specify initial concentrations for reactants, products, and a catalyst (case foe Ea analysis requires for very small concentration of a catalyst)
    c0_1 = {"CO": 1, "H2": 1, "ete": 1, "prod": 0, "I1_0L": 0.000001, "PMe3": 0.0000005}

    # Specify a reactant to study
    reactant_to_study = "PMe3"

    # Specify simulation time for apparent Ea analysis (in seconds)
    total_simulation_time_Ea = 10_000

    # Specify time in h at which rate of reactions (fluxes) is considered
    time = 2

    # Specify time step for the microkinetics simulations
    time_step = 1

    # Perform MicroKinetics Analysis
    analysis1 = ApparentEaAnalysis(temperature_value, 
                                   T_values_array_Ea, 
                                   reactant_concentration_array, 
                                   reactant_to_study,
                                   c0_1, 
                                   total_simulation_time_Ea, 
                                   time, 
                                   reactions, 
                                   time_step
    )

    # Specify reactant initial concentraiton from the concentration array that will be considered for ln(ri) vs. 1/T plot
    reactant_initial_concentration_plot = reactant_concentration_array[0]

    # Specify reactions that will be on plot of Ea vs. c0 (flux based)
    reactions_to_plot = [reaction for reaction in reactions if "prod" in reaction]

    # Specify compounds that will be on plot of Ea vs. c0 (compound based)
    compounds_to_plot = ["prod"]

    # Calculates and saves as csv dataframe with ri (flux) related parameters if wasn't calculated before for the specified parameters
    df_flux = analysis1.df_flux

    # Calculates and saves as csv dataframe with vi (rate) related parameters if wasn't calculated before for the specified parameters
    df_rate = analysis1.df_rate

    # Adjust nrows and ncols to your case
    analysis1.plot_ln_ri_vs_1_over_T(df_flux = df_flux, reactant_initial_concentration_plot = reactant_initial_concentration_plot, nrows = 5, ncols = 5, figsize = (15, 15))
    analysis1.plot_Ea_vs_c0_flux_based(df_flux = df_flux, reactions_to_plot = reactions_to_plot, nrows = 2, ncols = 3, figsize = (15, 10), log_x = True)
    analysis1.plot_Ea_vs_c0_rate_based(df_rate = df_rate, compounds_to_plot = compounds_to_plot, nrows = 1, ncols = 1, figsize = (15, 10), log_x = True)

    # Indicate the energy shift for the calculation of degree rate control (DRC)
    e_shift = 0.01

    # Specify range of T for which simulations will be calculated (Ea analysis)
    T_values_array_drc = np.linspace(start = temperature_value - 25,
                                     stop = temperature_value + 25,
                                     num = 1,
                                     endpoint = True
    )

    # Specify simulation time for apparent Ea analysis (in seconds)
    total_simulation_time_drc = 10_000

    # Specify number of cores to compute drc
    cores = 8

    analysis2 = DRCAnalysis(reactant_concentration_array,
                            T_values_array_drc,
                            total_simulation_time_drc,
                            c0_1,
                            reactant_to_study,
                            reactions,
                            e_shift, 
                            cores
    )

    # Temperature(s) that will be considered for drc vs. c0 plot (this variable can be an array if multiple plots are needed)
    temperature_value_drc_plot = T_values_array_drc[0]

    # Calculates and saves as csv dataframe with degree rate control (drc) related parameters if wasn't calculated before for the specified parameters
    df_drc = analysis2.df_drc

    # Adjust nrows and ncols to your case
    analysis2.plot_drc_vs_c0(nrows = 5, ncols = 5, figsize = (15, 15), df_drc = df_drc, temperature_value = temperature_value_drc_plot, log_x = True)

    # You can set now higher concentraiton of catalyst if needed
    c0_2 = {"CO": 0.1, "H2": 0.1, "ete": 0.1, "prod": 0, "I1_0L": 0.000005, "PMe3": 1}

    # Specify cycle names
    cycles = ["0L", "1L"]

    # Specify simulation time for "normal" microkinetics (in seconds)
    total_simulation_time_microkinetics = 100_000

    # Define reactants that could affect final concentration of a product (can be limiting reactant)
    main_reactants = ["CO", "H2", "ete"]

    # Specify compounds that will appear on a concentration evolution plot
    compounds_to_plot = ["I7_0L", "I7_1L", "I1_0L", "I1_1L"]
    
    # Set the time (in h) at which concentration of catalyst will be calculated (for C(catalyst) in each cycle in a system vs. initial concentration of studied reactant)
    catalyst_concentration_time = 1

    # Set the percentage of product conversion
    percentage_of_convertion = 0.99

    # Specify temperature(s) for concentration evolution plot from the temperature values array 
    temperatures = T_values_array_Ea[0]

    analysis3 = MicroKinetics(temperature_value, 
                              total_simulation_time_microkinetics, 
                              reactant_concentration_array, 
                              reactant_to_study, 
                              c0_2,
                              cycles, 
                              main_reactants, 
                              compounds_to_plot, 
                              catalyst_concentration_time, 
                              percentage_of_convertion,
                              time_step
    )

    # Adjust nrows and ncols to your case
    analysis3.plot_catalyst_concentration_vs_reactant(log_x = True, figsize = (15, 10))
    analysis3.plot_concentration_evolution(temperatures, log_y = True, log_x = True, nrows = 3, ncols = 3, figsize = (15, 15))
    analysis3.plot_reactant_vs_product_conversion(log_x = True, figsize = (15, 10))

    plt.show()

if __name__ == "__main__":
    main()