"""Barriers of each step from the Gibbs energies of its species (used by thermochemistry.py)"""

DIFFUSION_BARRIER = 4
G_COMPOUNDS_OUTPUT_DIR_NAME = "G_values_of_compounds"
REACTION_DF_OUTPUT_DIR_NAME = "G_values_of_reactions"


class ReactionFileParser:
    """Splits reaction strings into reactants and products"""

    @staticmethod
    def parse_reaction(reaction):
        """Parses a reaction string and returns the reactants and products"""
        reactants_str, products_str = reaction.split("=")
        reactants = [r.strip() for r in reactants_str.split("+")]
        products = [p.strip() for p in products_str.split("+")]
        return reactants, products


class GibbsEnergyCalculator:
    """Computes Gibbs energies (kcal/mol) of transition states and of direct/inverse reactions"""

    def __init__(self, compound_energy_dataframe):
        self.compound_energy_dataframe = compound_energy_dataframe

    def calculate_gibbs_free_energy_of_transition_state(self, ts, reactants, products):
        """Calculates Gibbs Free Energy of a transition state"""
        if ts == "-":
            return (
                max(
                    sum(
                        self.compound_energy_dataframe.loc[
                            reactant, "Gibbs Free Energies"
                        ]
                        for reactant in reactants
                    ),
                    sum(
                        self.compound_energy_dataframe.loc[
                            product, "Gibbs Free Energies"
                        ]
                        for product in products
                    ),
                )
                + DIFFUSION_BARRIER
            )
        else:
            return self.compound_energy_dataframe.loc[ts, "Gibbs Free Energies"]

    def calculate_gibbs_energy_of_reaction(self, compounds, ts_gibbs_free_energy):
        """Calculates Gibbs free energy for a direct/inverse reaction"""
        return ts_gibbs_free_energy - sum(
            self.compound_energy_dataframe.loc[compound, "Gibbs Free Energies"]
            for compound in compounds
        )

    def calculate_direct_inverse_reactions_gibbs_free_energies(self, reactions_df):
        """Calculates direct and inverse Gibbs Free Energies for each reaction"""
        Gdir_list = []
        Ginv_list = []

        for _, row in reactions_df.iterrows():
            reaction = row["Rx"]
            ts = row["TS"]
            # Get reactants and products of a reaction
            reactants, products = ReactionFileParser.parse_reaction(reaction)
            ts_gibbs_free_energy = self.calculate_gibbs_free_energy_of_transition_state(
                ts, reactants, products
            )

            # Calculate G for direct and inverse reactions
            Gdir = self.calculate_gibbs_energy_of_reaction(
                reactants, ts_gibbs_free_energy
            )
            Ginv = self.calculate_gibbs_energy_of_reaction(
                products, ts_gibbs_free_energy
            )

            Gdir_list.append(Gdir)
            Ginv_list.append(Ginv)

        return Gdir_list, Ginv_list
