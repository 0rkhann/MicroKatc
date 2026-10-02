"""Helper functions shared by the analyses (cache names, catalyst concentrations, conversion times, pressure)"""

import hashlib
import json

R_L_atm_per_mol_K = 0.082057366080960


class AuxiliaryFunctions:
    """Collection of static helper functions used across the analyses"""

    @staticmethod
    def inputs_hash(*inputs):
        """Short stable hash of the inputs, appended to cache file names so any changed input forces a recompute"""
        dumped = json.dumps(
            inputs,
            sort_keys=True,
            default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o),
        )
        return hashlib.md5(dumped.encode()).hexdigest()[:8]

    @staticmethod
    def get_concentration_for_cycle(simulations_dfs, time, intermediates):
        """Returns concentrations of a catalyst in a cycle at the specified time"""
        concentrations = []

        for df in simulations_dfs:
            # Find concentration of a catalyst in a cycle by summing concentrations of intermediates in this cycle
            concentration_sum = (
                df.loc[df["time"] == time, intermediates].sum(axis=1).item()
            )
            concentrations.append(concentration_sum)

        return concentrations

    @staticmethod
    def get_concentrations_of_catalyst(
        simulations_dfs, time, cycles_intermediates_dict
    ):
        """Returns concentrations of a catalyst in each cycle of a system at the specified time"""
        # Calculate catalyst concentrations for each cycle
        catalyst_concentrations_per_cycle = [
            AuxiliaryFunctions.get_concentration_for_cycle(
                simulations_dfs, time, intermediates
            )
            for intermediates in cycles_intermediates_dict.values()
        ]

        return catalyst_concentrations_per_cycle

    @staticmethod
    def compute_time_of_product_conversion_given_reactant_concentration(
        simulation_dfs,
        reactant_concentration_array,
        reactant_to_study,
        thresholds,
        product="prod",
    ):
        """First time each simulation's product concentration exceeds its threshold (one per simulation)"""
        times_conv_reac_conc = []

        for i, simulation_df in enumerate(simulation_dfs):
            try:
                # Find the first time value when concentration of a product is greater than a set product concentration threshold
                product_conversion_threshold_concentration = thresholds[i]
                filtered_df = simulation_df[
                    simulation_df[product] > product_conversion_threshold_concentration
                ]

                # First time at when product reach the product conversion threshold concentration
                first_convergence_time = filtered_df.iloc[0]["time"]

                # Access the 'time' column of the first row meeting the condition
                times_conv_reac_conc.append(
                    (first_convergence_time, reactant_concentration_array[i])
                )

            except Exception:
                print(
                    f"C({reactant_to_study}) = {reactant_concentration_array[i]}M couldn't reach the threshold of {product_conversion_threshold_concentration}M in the specified simulation time. It won't be included to the graph"
                )
                continue

        return times_conv_reac_conc

    @staticmethod
    def compute_pressure_value(temperature_value):
        """Computes pressure to satisfy liquid medium at standard conditions (1M)"""
        return 1 * R_L_atm_per_mol_K * temperature_value
