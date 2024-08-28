import pandas as pd
import re
import os
from file_operations import CURRENT_DIRECTORY
from calculating_G_for_microkinetics import G_COMPOUNDS_OUTPUT_DIR_NAME
import subprocess

R_L_atm_per_mol_K = 0.082057366080960

class AuxiliaryFunctions:
    @staticmethod
    def find_intermediates_of_cycle(cycles):
        """Returns intermediates belonging to a corresponding cycle as a dictionary"""
        compounds = pd.read_csv("compounds.csv")["Compounds"].to_list()
        cycles_intermediates_dict = {}

        for cycle in cycles:
            # Regex pattern for intermediates (I1_0L, I5_1L...)
            pattern = re.compile(f"I[\d]+_{cycle}")

            # Find intermediates matching the pattern for the current cycle
            intermediates = [intermediate for intermediate in compounds if pattern.match(intermediate)]
            cycles_intermediates_dict[cycle] = intermediates

        return cycles_intermediates_dict

    @staticmethod
    def find_limiting_reactant_concentration(initial_concentrations, reactants):
        """Returns concentration of a limiting reactant in a system"""
        min_initial_reactant_concentrations = min([initial_concentrations[reactant] for reactant in reactants])
        return min_initial_reactant_concentrations

    @staticmethod
    def get_concentration_for_cycle(simulations_dfs, time, intermediates):
        """Returns concentrations of a catalyst in a cycle at the specified time"""
        concentrations = []

        for df in simulations_dfs:
            # Find concentration of a catalyst in a cycle by summing concentrations of intermediates in this cycle
            concentration_sum = df.loc[df['time'] == time, intermediates].sum(axis = 1).item()
            concentrations.append(concentration_sum)

        return concentrations

    @staticmethod
    def get_concentrations_of_catalyst(simulations_dfs, time, cycles_intermediates_dict):
        """Returns concentrations of a catalyst in each cycle of a system at the specified time"""
        # Calculate catalyst concentrations for each cycle
        catalyst_concentrations_per_cycle = [
            AuxiliaryFunctions.get_concentration_for_cycle(simulations_dfs, time, intermediates)
            for intermediates in cycles_intermediates_dict.values()
        ]

        return catalyst_concentrations_per_cycle

    @staticmethod
    def compute_time_of_product_conversion_given_reactant_concentration(simulation_dfs, reactant_concentration_array, reactant_to_study, product_conversion_threshold_concentration):
        """Returns the first times when concentration of a product is greater than a product concentration threshold"""
        times_conv_reac_conc = []

        for i, simulation_df in enumerate(simulation_dfs):
            
            try:
                # Find the first time value when concentration of a product is greater than a set product concentration threshold
                filtered_df = simulation_df.query("prod > @product_conversion_threshold_concentration")

                # First time at when product reach the product conversion threshold concentration
                first_convergence_time = filtered_df.iloc[0]['time']

                # Access the 'time' column of the first row meeting the condition
                times_conv_reac_conc.append((first_convergence_time, reactant_concentration_array[i]))
            
            except:
                print(f"C({reactant_to_study}) = {reactant_concentration_array[i]}M couldn't reach the threshold of {product_conversion_threshold_concentration}M in the specified simulation time. It won't be included to the graph")
                continue

        return times_conv_reac_conc
    
    @staticmethod
    def reactions_number(reaction_df):
        """Returns reactions as a list and number of reactions"""
        reactions = reaction_df["Rx"].to_list()
        return reactions

    @staticmethod
    def calculate_G_values(temperature, pressure):
        """Calculates G of compounds and reactions in a system at specified T and P"""
        reaction_df_file_name = f"reaction_df_{temperature}K_{pressure:.5e}atm.csv"
        G_compounds_file_name = f"G_values_at_{temperature}K_{pressure:.5e}atm.csv"

        G_compounds_file_path = os.path.join(CURRENT_DIRECTORY, G_COMPOUNDS_OUTPUT_DIR_NAME, G_compounds_file_name)

        # Calculations will be executed if the file with specified T and P doesn't in a dedicated directory
        if not os.path.exists(G_compounds_file_path): 
            # Path to the Bash script
            script_path = os.path.join(CURRENT_DIRECTORY, "get_G_compounds.sh")
            print(f"Creating and saving: {reaction_df_file_name}")
            print(f"Creating and saving: {G_compounds_file_name}")
            subprocess.run(['bash', script_path, f"{temperature}", f"{pressure}"], capture_output = True, text = True)

    @staticmethod
    def compute_pressure_value(temperature_value):
        """Computes pressure to satisfy liquid medium at standard conditions (1M)"""
        return 1 * R_L_atm_per_mol_K * temperature_value