import os
import pandas as pd
from file_operations import FileOperations, CURRENT_DIRECTORY
from bash_parsing import BashParser

DIFFUSION_BARRIER = 4 
G_COMPOUNDS_OUTPUT_DIR_NAME = "G_values_of_compounds"
REACTION_DF_OUTPUT_DIR_NAME = "G_values_of_reactions"

class TemperaturePressureParser:
    @staticmethod
    def get_temperature_pressure():
        """Fetch temperature and pressure values from bash input"""
        return BashParser.get_temperature_pressure_from_input()

class ReactionFileParser:
    def __init__(self, reaction_file="reactions.csv"):
        self.reactions_df = pd.read_csv(reaction_file, sep = ",")

    @staticmethod
    def parse_reaction(reaction):
        """Parses a reaction string and returns the reactants and products"""
        reactants_str, products_str = reaction.split('=')
        reactants = [r.strip() for r in reactants_str.split('+')]
        products = [p.strip() for p in products_str.split('+')]
        return reactants, products

class GibbsEnergyCalculator:
    def __init__(self, compound_energy_dataframe):
        self.compound_energy_dataframe = compound_energy_dataframe

    def calculate_gibbs_free_energy_of_transition_state(self, ts, reactants, products):
        """Calculates Gibbs Free Energy of a transition state"""
        if ts == "-":
            return max(sum(self.compound_energy_dataframe.loc[reactant, "Gibbs Free Energies"] for reactant in reactants),
                       sum(self.compound_energy_dataframe.loc[product, "Gibbs Free Energies"] for product in products)) + DIFFUSION_BARRIER
        else:
            return self.compound_energy_dataframe.loc[ts, "Gibbs Free Energies"]

    def calculate_gibbs_energy_of_reaction(self, compounds, ts_gibbs_free_energy):
        """Calculates Gibbs free energy for a direct/inverse reaction"""
        return ts_gibbs_free_energy - sum(self.compound_energy_dataframe.loc[compound, "Gibbs Free Energies"] for compound in compounds)

    def calculate_direct_inverse_reactions_gibbs_free_energies(self, reactions_df):
        """Calculates direct and inverse Gibbs Free Energies for each reaction"""
        Gdir_list = []
        Ginv_list = []

        for _, row in reactions_df.iterrows():
            reaction = row['Rx']
            ts = row['TS']
            # Get reactants and products of a reaction
            reactants, products = ReactionFileParser.parse_reaction(reaction)
            ts_gibbs_free_energy = self.calculate_gibbs_free_energy_of_transition_state(ts, reactants, products)

            # Calculate G for direct and inverse reactions
            Gdir = self.calculate_gibbs_energy_of_reaction(reactants, ts_gibbs_free_energy)
            Ginv = self.calculate_gibbs_energy_of_reaction(products, ts_gibbs_free_energy)

            Gdir_list.append(Gdir)
            Ginv_list.append(Ginv)

        return Gdir_list, Ginv_list

class CompoundsDataHandler:
    def __init__(self, df_compounds_file_name="compounds.csv"):
        self.df_compounds_file_name = df_compounds_file_name

    def create_compounds_df(self, compounds):
        """Creates and saves dataframe with compounds in the system if it doesn't exist"""
        compounds_df_path = os.path.join(CURRENT_DIRECTORY, self.df_compounds_file_name)
        if not os.path.exists(compounds_df_path):
            print("Creating compounds.csv file...")
            compounds.to_csv(f"{self.df_compounds_file_name}", index=False)

    def save_compounds_G_values(self, compound_energy_dataframe, G_compounds_file_name):
        """Saves compounds data to CSV and move to the appropriate directory"""
        os.makedirs(G_COMPOUNDS_OUTPUT_DIR_NAME, exist_ok=True)
        compound_energy_dataframe.to_csv(G_compounds_file_name, index=True)
        FileOperations.move_to_output_directory(G_COMPOUNDS_OUTPUT_DIR_NAME, G_compounds_file_name)

class ReactionDataHandler:
    def save_reaction_data(self, reactions_df, reaction_df_file_name):
        """Saves reaction data to CSV and move to the appropriate directory"""
        os.makedirs(REACTION_DF_OUTPUT_DIR_NAME, exist_ok=True)
        reactions_df.to_csv(reaction_df_file_name, index=False)
        FileOperations.move_to_output_directory(REACTION_DF_OUTPUT_DIR_NAME, reaction_df_file_name)

class ReactionCalculator:
    def __init__(self, reaction_file="reactions.csv"):
        # Initialize parsers and handlers
        self.temperature, self.pressure = TemperaturePressureParser.get_temperature_pressure()
        self.reaction_file_parser = ReactionFileParser(reaction_file)
        self.compounds_data_handler = CompoundsDataHandler()
        self.reaction_data_handler = ReactionDataHandler()

        self.G_compounds_file_name = f"G_values_at_{self.temperature}K_{self.pressure:.5e}atm.csv"
        self.reaction_df_file_name = f"reaction_df_{self.temperature}K_{self.pressure:.5e}atm.csv"

        # Fetch compound data and calculate Gibbs free energies
        self.compound_energy_dataframe, self.compounds = BashParser.parse_output_from_bash_script()
        self.compounds_data_handler.create_compounds_df(self.compounds)
        self.compounds_data_handler.save_compounds_G_values(self.compound_energy_dataframe, self.G_compounds_file_name)

        # Calculate Gibbs free energy for reactions
        self.gibbs_calculator = GibbsEnergyCalculator(self.compound_energy_dataframe)
        Gdir_values, Ginv_values = self.gibbs_calculator.calculate_direct_inverse_reactions_gibbs_free_energies(self.reaction_file_parser.reactions_df)

        # Update reaction dataframe and save
        self.reaction_file_parser.reactions_df['Gdir'], self.reaction_file_parser.reactions_df['Ginv'] = Gdir_values, Ginv_values
        self.reaction_data_handler.save_reaction_data(self.reaction_file_parser.reactions_df, self.reaction_df_file_name)

if __name__ == "__main__":
    ReactionCalculator()