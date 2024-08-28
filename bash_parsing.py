import sys
import pandas as pd

HARTREE_TO_KCAL_MOL = 627.509

class BashParser:
    @staticmethod
    def get_temperature_pressure_from_input():
        """Get temperature and pressure from sys.argv"""
        temperature = float(sys.argv[2])
        pressure = float(sys.argv[3])
        return temperature, pressure

    @staticmethod
    def parse_output_from_bash_script():
        """Parses output from a bash script and creates a dataframe with compounds and corresponding to them G and creates dataframe with compounds"""
        df = pd.read_csv("energies.temp", sep = " ", header = None)
        df = df.rename(columns = {0: "Compounds", 1: "Gibbs Free Energies"})
        df_compounds = df["Compounds"]
        df = df.set_index("Compounds")
        df['Gibbs Free Energies'] = df['Gibbs Free Energies'].astype(float) * HARTREE_TO_KCAL_MOL
        return df, df_compounds