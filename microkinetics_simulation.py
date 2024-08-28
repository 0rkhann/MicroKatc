import copasi_parser as cpx
from plotting_functions import PlotFunctions
from auxilary_functions import AuxiliaryFunctions
from file_operations import FileOperations, CURRENT_DIRECTORY
from calculating_G_for_microkinetics import REACTION_DF_OUTPUT_DIR_NAME
import os
from collections.abc import Iterable

CONVERT_SECONDS_TO_HOURS = 1 / 3600
SIMULATIONS_OUTPUT_DIR_NAME = "microkinetics_simulations"

class SimulationHandler:
    def __init__(self, T, P, reactant_to_study, total_simulation_time, time_step = 1):
        self.temperature_value = T
        self.pressure_value = P
        self.reactant_to_study = reactant_to_study
        self.total_simulation_time = total_simulation_time
        self.time_step = time_step
        self.datafile = os.path.join(CURRENT_DIRECTORY, REACTION_DF_OUTPUT_DIR_NAME, f"reaction_df_{self.temperature_value}K_{self.pressure_value:.5e}atm.csv")

    def run_simulation(self, simulation_filename, c0, convergence_flag=True):
        """Runs the microkinetics simulation and returns the results as a DataFrame"""
        ch = cpx.prepare_copasi_model(reactions_file = self.datafile, temp = self.temperature_value, initial_concentrations = c0, csv_delim = ",")
        traj_task, result_flag = cpx.time_course_simulation(ch, total_time = self.total_simulation_time, time_step = self.time_step,
                                                            delay_run = False, model_file = "model.cps", outfile = simulation_filename, 
                                                            add_fluxes = True, add_rates = True)
        simulation_df = cpx.read_simulation(simulation_filename)

        if convergence_flag:
            simulation_df = self.check_copasi_convergence(simulation_df, simulation_filename, c0)

        simulation_df.loc[:, "time"] *= CONVERT_SECONDS_TO_HOURS
        ch.self_destruct()
        return simulation_df

    def check_copasi_convergence(self, simulation_df, simulation_filename, c0):
        """COPASI rarely may not converge to a solution due to complexity of differential equations nature. It checks whether a simulation is converged"""
        convergence_counter = 0
        while simulation_df["time"].iloc[-1] != self.total_simulation_time:
            print(f"{simulation_filename} couldn't converge to a solution. Resimulating...")
            simulation_df = self.run_simulation(simulation_filename, c0, False)
            FileOperations.move_to_output_directory(SIMULATIONS_OUTPUT_DIR_NAME, simulation_filename)
            convergence_counter += 1
            if convergence_counter == 5:
                raise ConvergenceError(f"Solution can't be found for {simulation_filename}. Please, change the input conditions")
                
        return simulation_df
            
    def get_simulation_df(self, c0):
        """Checks if the simulation results exist and returns them, or runs a new simulation if needed"""
        os.makedirs(SIMULATIONS_OUTPUT_DIR_NAME, exist_ok = True)
        simulation_filename = "_".join([
                                       f"{reactant}_{concentration:.6e}M" if reactant == self.reactant_to_study else f"{reactant}_{concentration}M"
                                       for reactant, concentration in c0.items()
                                      ]) \
                                        + f"_{self.temperature_value}K_{self.pressure_value:.5e}atm" + f"_{self.total_simulation_time}s"
        
        simulation_file_path = os.path.join(CURRENT_DIRECTORY, SIMULATIONS_OUTPUT_DIR_NAME, simulation_filename)

        if os.path.exists(simulation_file_path):
            print(f"Reading: {simulation_filename}")
            simulation_df = cpx.read_simulation(simulation_file_path)
            simulation_df.loc[:, "time"] *= CONVERT_SECONDS_TO_HOURS
            return simulation_df
        else:
            print(f"Simulating: {simulation_filename}")
            simulation_df = self.run_simulation(simulation_filename, c0)
            FileOperations.move_to_output_directory(SIMULATIONS_OUTPUT_DIR_NAME, simulation_filename)
            return simulation_df

class ConvergenceError(Exception):
    pass

class MicroKinetics:
    def __init__(self, temperature_value, total_simulation_time, reactant_concentration_array, reactant_to_study, 
                 c0, cycles, main_reactants, compounds_to_plot, catalyst_concentration_time, percentage_of_convertion, time_step):
        
        self.temperature_value = temperature_value
        self.pressure_value = AuxiliaryFunctions.compute_pressure_value(self.temperature_value)
        self.reactant_to_study = reactant_to_study
        self.total_simulation_time = total_simulation_time
        self.c0 = c0
        self.cycles = cycles
        self.main_reactants = main_reactants
        self.compounds_to_plot = compounds_to_plot
        self.reactant_concentration_array = reactant_concentration_array
        self.catalyst_concentration_time = catalyst_concentration_time
        self.percentage_of_convertion = percentage_of_convertion
        
        self.simulation_handler = SimulationHandler(self.temperature_value, 
                                                    self.pressure_value, 
                                                    self.reactant_to_study, 
                                                    self.total_simulation_time, 
                                                    time_step = time_step
        )
        self._simulations_dfs = None

        self.plot_functions = PlotFunctions()

    @property
    def simulations_dfs(self):
        if self._simulations_dfs is None:
            print("Starting computing simulation dataframes...")
            # List of simulation dataframes for different initial concentrations of reactant_to_study
            self._simulations_dfs = [
                self.simulation_handler.get_simulation_df({**self.c0, self.reactant_to_study: initial_concentration})
                for initial_concentration in self.reactant_concentration_array
            ]
        return self._simulations_dfs

    def plot_concentration_evolution(self, temperatures, log_y, log_x, nrows, ncols, figsize):
        if not isinstance(temperatures, Iterable):
            temperatures = [temperatures] 

        for temperature in temperatures:
            fig_name = f"C(compounds)_C({self.reactant_to_study})_{temperature}K_{self.pressure_value:.5e}atm_{self.compounds_to_plot}.svg"
            fig = self.plot_functions.plot_concentrations_evolution(
                                                self.reactant_to_study,
                                                self.simulations_dfs,
                                                self.reactant_concentration_array,
                                                self.compounds_to_plot,
                                                nrows,
                                                ncols,
                                                figsize,
                                                fig_name,
                                                temperature,
                                                log_y,
                                                log_x
            )
            self.plot_functions.save_figure(fig_name)

    def plot_catalyst_concentration_vs_reactant(self, log_x, figsize):
        fig_name = f"C(catalyst)_C({self.reactant_to_study})_{self.temperature_value}K_{self.pressure_value:.5e}atm.svg"
        cycles_intermediates_dict = AuxiliaryFunctions.find_intermediates_of_cycle(self.cycles)
        catalyst_concentrations_per_cycle = AuxiliaryFunctions.get_concentrations_of_catalyst(self.simulations_dfs, self.catalyst_concentration_time, cycles_intermediates_dict)

        fig = self.plot_functions.plot_concentration_of_catalyst_versus_studied_reactant(self.reactant_concentration_array, 
                                                                                         catalyst_concentrations_per_cycle, 
                                                                                         self.cycles, 
                                                                                         self.reactant_to_study, 
                                                                                         figsize, 
                                                                                         fig_name,
                                                                                         log_x
        )
        self.plot_functions.save_figure(fig_name)

    def plot_reactant_vs_product_conversion(self, log_x, figsize):
        fig_name = f"C({self.reactant_to_study})_total_t_product_conversion_{self.temperature_value}K_{self.pressure_value:.5e}atm.svg"
        product_conversion_threshold_concentration = AuxiliaryFunctions.find_limiting_reactant_concentration(self.c0, self.main_reactants) * self.percentage_of_convertion

        times_conv_reac_conc = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(self.simulations_dfs, 
                                                                                                                  self.reactant_concentration_array, 
                                                                                                                  self.reactant_to_study, 
                                                                                                                  product_conversion_threshold_concentration)

        fig = self.plot_functions.plot_reactant_concentration_vs_product_conversion(times_conv_reac_conc,
                                                                                    self.reactant_to_study, 
                                                                                    figsize, 
                                                                                    fig_name,
                                                                                    log_x
        )
        self.plot_functions.save_figure(fig_name)