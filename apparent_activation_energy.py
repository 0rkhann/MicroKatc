from microkinetics_simulation import SimulationHandler, SIMULATIONS_OUTPUT_DIR_NAME
import numpy as np
import pandas as pd
from auxilary_functions import AuxiliaryFunctions
from plotting_functions import PlotFunctions
from scipy.constants import R
import copasi_parser as cpx
import os
from file_operations import CURRENT_DIRECTORY, FileOperations
from calculating_G_for_microkinetics import REACTION_DF_OUTPUT_DIR_NAME

J_TO_KCAL = 4184

class ReactionDataHandler:
    """Handles the retrieval and storage of reaction data"""

    @staticmethod
    def get_dataframe_with_data(T_values_array, time, reactant_to_study, reactant_concentration_array, c0, total_simulation_time, time_step):
        """Calculates ln(ri) and ln(vi) values at a specified time for each initial concentration of a studied reactant and returns a DataFrame"""
        data = []
        fluxes, rates = [], []

        # Iterate over each initial concentration and temperature
        for initial_concentration in reactant_concentration_array:
            c0_copy = c0.copy()
            c0_copy[reactant_to_study] = initial_concentration

            for T_value in T_values_array:
                pressure_value = AuxiliaryFunctions.compute_pressure_value(T_value)
                AuxiliaryFunctions.calculate_G_values(T_value, pressure_value)

                simulation_calculator = SimulationHandler(T_value, pressure_value, reactant_to_study, total_simulation_time, time_step)
                simulation_df = simulation_calculator.get_simulation_df(c0_copy)

                # Do operation only once as columns are the same each iteration
                if fluxes == []:
                    fluxes = [flux_name for flux_name in simulation_df if flux_name.endswith(".Flux")]
                    rates = [rate_name for rate_name in simulation_df if rate_name.endswith(".Rate")]

                # Collect flux data
                data.extend(ReactionDataHandler.collect_log_values(simulation_df, 
                                                                   fluxes, 
                                                                   time, 
                                                                   initial_concentration, 
                                                                   T_value, 
                                                                   reactant_to_study)
                )

                # Collect rate data
                data.extend(ReactionDataHandler.collect_log_values(simulation_df, 
                                                                   rates, 
                                                                   time, 
                                                                   initial_concentration, 
                                                                   T_value, 
                                                                   reactant_to_study)
                )

        return pd.DataFrame(data)
    
    @staticmethod
    def collect_log_values(simulation_df, keys, time, initial_concentration, T_value, reactant_to_study):
        """Collects and logs values (ri or vi) for specified keys in the simulation data"""
        collected_data = []
        for key in keys:
            value = simulation_df.loc[simulation_df['time'] == time, key].item()
            log_value = np.log(abs(value))

            data_entry = {
                "name": key,
                f"c0({reactant_to_study})": initial_concentration,
                "T": T_value,
                "1/T": 1 / T_value,
                f"ln_value": log_value
            }
            collected_data.append(data_entry)

        return collected_data

class ReactionParameterCalculator:
    """Calculates reaction parameters like activation energy and R2"""

    @staticmethod
    def calculate_reaction_parameters(df, T_range, reactant_concentration_array, reactant_to_study, calculation_type, reactions = None):
        """Calculates apparent activation energy from the slope: ln(ri) = -Ea/RT + ln(A) and ln(vi) = -Ea/RT + ln(A). Returns a dataframe"""
        data = []
        keys = [key for key in df["name"].unique() if key.endswith(".Flux")] if calculation_type == "ri" else [key for key in df["name"].unique() if key.endswith(".Rate")]

        for reactant_c0 in reactant_concentration_array:
            for i, key in enumerate(keys):
                # Dataframe with data of current flux (reaction) and current initial concentration 
                filtered_df = df[(df["name"] == key) & (np.isclose(df[f"c0({reactant_to_study})"], reactant_c0, atol = 0.0))]
                
                if not filtered_df.empty:
                    x_values = filtered_df["1/T"].values
                    y_values = filtered_df["ln_value"].values
                    
                    slope, intercept = np.polyfit(x_values, y_values, 1)
                    
                    y_pred = slope * x_values + intercept
                    ss_res = np.sum((y_values - y_pred) ** 2)
                    ss_tot = np.sum((y_values - np.mean(y_values)) ** 2)
                    r2 = 1 - (ss_res / ss_tot)

                    Ea = -slope * R / J_TO_KCAL
                    
                    data.append({
                            "name": key,
                            f"c0({reactant_to_study})": reactant_c0,
                            "R2": r2,
                            "Ea": Ea,
                            "slope": slope,
                            "intercept": intercept
                    } | {f"1/T({i})": x_values[i] for i in range(len(x_values))}
                      | {f"ln_value({i})": y_values[i] for i in range(len(y_values))}
                      | {f"T({i})": T_range[i] for i in range(len(T_range))}
                    )

                    if calculation_type == "ri":
                        data[-1]["reaction"] = reactions[i]
                    else:
                        data[-1]["compound"] = key.replace(".Rate", "")

        return pd.DataFrame(data)

class ApparentEaAnalysis:
    """Orchestrates the microkinetic analysis process including data handling, parameter calculations, and plotting"""

    def __init__(self, temperature_value, T_values_array, reactant_concentration_array, reactant_to_study, c0, total_simulation_time,
                 time, reactions, time_step):
        self.temperature_value = temperature_value
        self.T_values_array = T_values_array
        self.reactant_concentration_array = reactant_concentration_array
        self.reactant_to_study = reactant_to_study
        self.c0 = c0
        self.total_simulation_time = total_simulation_time
        self.time = time
        self.reactions = reactions
        self.time_step = time_step

        self.plot_function = PlotFunctions()

        self.df_flux_filename = f"df_flux_T_range_{self.T_values_array[0]}K_{self.T_values_array[-1]}K_C({self.reactant_to_study})_{self.reactant_concentration_array[0]}M_{self.reactant_concentration_array[-1]}M.csv"
        self.df_rate_filename = f"df_rate_T_range_{self.T_values_array[0]}K_{self.T_values_array[-1]}K_C({self.reactant_to_study})_{self.reactant_concentration_array[0]}M_{self.reactant_concentration_array[-1]}M.csv"

        self._df = None

        try:
            df_flux_path = os.path.join(CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, self.df_flux_filename)
            self._df_flux = pd.read_csv(df_flux_path)
            print(f"Reading existing flux dataframe with provide T range {self.T_values_array[::len(self.T_values_array)-1]} and reactant "
                  f"({self.reactant_to_study}) concentration range {self.reactant_concentration_array[::len(self.reactant_concentration_array)-1]}...")
        except:
            self._df_flux = None
        
        try:
            df_rate_path = os.path.join(CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, self.df_rate_filename)
            self._df_rate = pd.read_csv(df_rate_path)
            print(f"Reading existing rate dataframe with provide T range {self.T_values_array[::len(self.T_values_array)-1]} and reactant "
                  f"({self.reactant_to_study}) concentration range {self.reactant_concentration_array[::len(self.reactant_concentration_array)-1]}...")
        except:
            self._df_rate = None

    @property
    def df(self):
        if self._df is None:
            print("Computing dataframe with ln(ri) and ln(vi) values...")
            self._df = ReactionDataHandler.get_dataframe_with_data(self.T_values_array, 
                                                                   self.time, 
                                                                   self.reactant_to_study, 
                                                                   self.reactant_concentration_array, 
                                                                   self.c0, 
                                                                   self.total_simulation_time,
                                                                   self.time_step
            )
        return self._df

    @property
    def df_flux(self):
        if self._df_flux is None:
            print("Computing flux-based parameters...")
            self._df_flux = ReactionParameterCalculator.calculate_reaction_parameters(self.df, 
                                                                                      self.T_values_array, 
                                                                                      self.reactant_concentration_array, 
                                                                                      self.reactant_to_study, 
                                                                                      calculation_type = "ri", 
                                                                                      reactions = self.reactions
            )
            print("Saving df flux csv...")
            self._df_flux.to_csv(self.df_flux_filename, index = False)
            FileOperations.move_to_output_directory(SIMULATIONS_OUTPUT_DIR_NAME, self.df_flux_filename)
        return self._df_flux

    @property
    def df_rate(self):
        if self._df_rate is None:
            print("Computing rate-based parameters...")
            self._df_rate = ReactionParameterCalculator.calculate_reaction_parameters(self.df, 
                                                                                      self.T_values_array, 
                                                                                      self.reactant_concentration_array, 
                                                                                      self.reactant_to_study, 
                                                                                      calculation_type = "vi"
            )
            print("Saving df rate csv...")
            self._df_rate.to_csv(self.df_rate_filename, index = False)
            FileOperations.move_to_output_directory(SIMULATIONS_OUTPUT_DIR_NAME, self.df_rate_filename)
        return self._df_rate

    def plot_ln_ri_vs_1_over_T(self, df_flux, reactant_initial_concentration_plot, nrows, ncols, figsize):
        fig_name = f"ln(ri)_1_T_{self.reactant_to_study}_{reactant_initial_concentration_plot:.6e}_T_range_{self.T_values_array[0]}K_{self.T_values_array[-1]}K.svg"

        self.plot_function.plot_ln_ri_vs_1_over_T(nrows, 
                                                  ncols, 
                                                  df_flux, 
                                                  self.T_values_array,
                                                  figsize, 
                                                  reactant_initial_concentration_plot, 
                                                  self.reactant_to_study, 
                                                  fig_name,
        )
        self.plot_function.save_figure(fig_name)

    def plot_Ea_vs_c0_flux_based(self, df_flux, reactions_to_plot, nrows, ncols, figsize, log_x):
        fig_name = f"Ea_c0({self.reactant_to_study})_flux_based.svg"
        self.plot_function.plot_Ea_vs_reactant_c0(nrows, 
                                                  ncols, 
                                                  df_flux, 
                                                  figsize, 
                                                  self.reactant_to_study, 
                                                  reactions_to_plot, 
                                                  fig_name, 
                                                  log_x = log_x,
                                                  value_type = "ri"
        )
        self.plot_function.save_figure(fig_name)

    def plot_Ea_vs_c0_rate_based(self, df_rate, compounds_to_plot, nrows, ncols, figsize, log_x):
        fig_name = f"Ea_c0({self.reactant_to_study})_rate_based.svg"
        self.plot_function.plot_Ea_vs_reactant_c0(nrows, 
                                                  ncols, 
                                                  df_rate, 
                                                  figsize, 
                                                  self.reactant_to_study, 
                                                  compounds_to_plot, 
                                                  fig_name, 
                                                  log_x = log_x,
                                                  value_type = "vi"
        )
        self.plot_function.save_figure(fig_name)

class DRCAnalysis:
    def __init__(self, reactant_concentration_array, T_values_array, total_simulation_time, c0, reactant_to_study, reactions, e_shift, cores=1):
        self.reactant_concentration_array = reactant_concentration_array
        self.T_values_array = T_values_array
        self.total_simulation_time = total_simulation_time
        self.c0 = c0
        self.reactant_to_study = reactant_to_study
        self.reactions = reactions
        self.e_shift = e_shift
        self.cores = cores

        self.plot_function = PlotFunctions()

        self.df_drc_filename = f"df_drc_T_range_{self.T_values_array[0]}K_{self.T_values_array[-1]}K_C({self.reactant_to_study})_{self.reactant_concentration_array[0]}M_{self.reactant_concentration_array[-1]}M.csv"

        try:
            df_drc_path = os.path.join(CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, self.df_drc_filename)
            self._df_drc = pd.read_csv(df_drc_path)
            print(f"Reading existing drc dataframe with provide T range {self.T_values_array[::len(self.T_values_array)-1]} and reactant "
                  f"({self.reactant_to_study}) concentration range {self.reactant_concentration_array[::len(self.reactant_concentration_array)-1]}...")
        except:
            self._df_drc = None

    @property
    def df_drc(self):
        if self._df_drc is None:
            self._df_drc = self.calculate_degree_of_rate_control()
            print("Saving df drc csv...")
            self._df_drc.to_csv(self.df_drc_filename, index = False)
            FileOperations.move_to_output_directory(SIMULATIONS_OUTPUT_DIR_NAME, self.df_drc_filename)
            
        return self._df_drc

    def calculate_degree_of_rate_control(self):
        """Calculates degree of rate control coefficients for each step in a system given different temperatures and different initial concentraitons of reactant to study"""
        data = []

        for T_value in self.T_values_array:
            pressure_value = AuxiliaryFunctions.compute_pressure_value(T_value)
            datafile = os.path.join(CURRENT_DIRECTORY, REACTION_DF_OUTPUT_DIR_NAME, f"reaction_df_{T_value}K_{pressure_value:.5e}atm.csv")

            for initial_concentration in self.reactant_concentration_array:
                c0_copy = self.c0.copy()
                c0_copy[self.reactant_to_study] = initial_concentration

                print(f"Computing degree of rate control for {T_value}K, {pressure_value:.2e}atm, c({self.reactant_to_study}) = {initial_concentration}M...")

                cpx_model = cpx.prepare_copasi_model(reactions_file = datafile, 
                                                    temp = T_value, 
                                                    initial_concentrations = c0_copy, 
                                                    csv_delim = ","
                )

                drc_coefficients, _ = cpx.drc_calc(base_model = cpx_model,
                                                    temp = T_value,
                                                    c0 = c0_copy, 
                                                    total_time = self.total_simulation_time,
                                                    time_step = 1,
                                                    target_time = self.total_simulation_time / 2,
                                                    target_spc = "prod",
                                                    e_shift = self.e_shift,
                                                    cores = self.cores
                                                    
                )
                
                data.append({
                        f"c0({self.reactant_to_study})": initial_concentration,
                        "temperature": T_value,
                        **{reaction: coeff for reaction, coeff in zip(self.reactions, drc_coefficients)}
                })

                cpx_model.self_destruct()

        return pd.DataFrame(data)

    def plot_drc_vs_c0(self, nrows, ncols, figsize, df_drc, temperature_value, log_x):
        fig_name = f"drc_c0_{self.reactant_to_study}_{self.reactant_concentration_array[0]}_{self.reactant_concentration_array[-1]}.svg"
        self.plot_function.plot_drc_vs_c0(nrows, 
                                          ncols, 
                                          df_drc, 
                                          figsize, 
                                          self.reactions, 
                                          self.reactant_to_study, 
                                          fig_name, 
                                          temperature_value,
                                          log_x = log_x
        )
        self.plot_function.save_figure(fig_name)