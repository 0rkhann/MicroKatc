import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 
import os
import numpy as np
from file_operations import FileOperations

IMAGES_SIMULATIONS_OUTPUT_DIR_NAME = "microkinetics_simulations_images"

class PlotFunctions:
    def __init__(self):
        os.makedirs(IMAGES_SIMULATIONS_OUTPUT_DIR_NAME, exist_ok = True)

    @staticmethod
    def hide_unused_subplots(fig, axes):
        """Hides any unused subplots (if any)"""
        for j in range(len(axes)):
            if not axes[j].has_data(): 
                fig.delaxes(axes[j])

    def plot_concentration_of_catalyst_versus_studied_reactant(self, reactant_concentrations, catalyst_concentrations_per_cycle, cycles, reactant_to_study, figsize, fig_name, log_x=False):      
        """Plots concentration of catalysts in all cycles at a specific time for different initial concentrations of reactant to see how initial concentration
        of a reactant affects the conversion of a catalyst from one form to another (advancement from one cycle to another)"""
        print(f"Plotting: {fig_name}...")
        fig, axes = self.initialize_fig_axes_objects(figsize = figsize)

        # Iterate though each cycle and plot on the same graph catalyst concentration in each in a function of initial reactant concentration (at specified time)
        for ax in axes:
            for i, cycle in enumerate(cycles):
                self.plot_axe_from_arrays(x = reactant_concentrations, 
                                          y = catalyst_concentrations_per_cycle[i],
                                          xlabel = f"C({reactant_to_study}) (M)", 
                                          ylabel = f"C(Catalyst) (M)", 
                                          title = f"C(catalyst) in different cycles vs. C({reactant_to_study})", 
                                          ax = ax,
                                          legend = f"Catalyst in {cycle} cycle"
                )

            if log_x:
                ax.set_xscale("log")

        self.hide_unused_subplots(fig, axes)
        return fig

    def plot_reactant_concentration_vs_product_conversion(self, times_conv_reac_conc, reactant_to_study, figsize, fig_name, log_x=False):
        """Plots graph of time of total product conversion in a function of reactant concentration to see how initial concentration
        of a reactant affects the time of getting 99% of theoretically possible quantity of a product"""
        print(f"Plotting: {fig_name}...")
        fig, axes = self.initialize_fig_axes_objects(figsize = figsize)

        times_of_product_conversion = [x_y_pair[0] for x_y_pair in times_conv_reac_conc]
        reactant_concentrations =  [x_y_pair[1] for x_y_pair in times_conv_reac_conc]

        for ax in axes:
            self.plot_axe_from_arrays(x = reactant_concentrations,
                                      y = times_of_product_conversion,
                                      xlabel = f"C({reactant_to_study}) (M)", 
                                      ylabel = f"Time of product conversion",
                                      title = f"Time of product conversion vs. C({reactant_to_study})",
                                      ax = ax
            )

            if log_x:
                ax.set_xscale("log")

        self.hide_unused_subplots(fig, axes)
        return fig

    @staticmethod
    def initialize_fig_axes_objects(nrows=1, ncols=1, figsize=(15, 15)):
        """Creates a fig and ax objects for plotting considering the number of rows and columns"""
        fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = (figsize[0] * ncols * 0.4, figsize[1] * nrows * 0.4))
        axes = axes.flatten() if nrows * ncols > 1 else [axes]
        fig.subplots_adjust(hspace = 0.4, wspace = 0.4) if ncols > 1 else None
        return fig, axes

    @staticmethod
    def plot_axe_from_arrays(x, y, xlabel, ylabel, title, ax, legend=None):
        """Plots a graph with provided x and y values. Plots with "o" markers on it and as a continious line"""
        ax.plot(x, y, "o", label = legend, linestyle = "-")
        if legend:
            ax.legend()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    @staticmethod
    def plot_axe_from_df_seaborn_groupplot(x, y, xlabel, ylabel, title, ax, hue=None, plot_df=None, legend=False, palette=None, log_y=False, log_x = False):
        """Plots a graph with provided x and y (used in case of complex data coming from dataframe)"""
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        if log_y:
            ax.set_yscale("log")

        if log_x:
            ax.set_xscale("log")

        sns.lineplot(data = plot_df,
                     x = x,
                     y = y,
                     hue = hue,
                     palette = palette,
                     legend = legend,
                     ax = ax,
                     lw = 3
        )

    def plot_concentrations_evolution(self, reactant_to_study, simulations_dfs, reactant_concentration_array, compounds_for_plotting, nrows, ncols, figsize, fig_name, temperature, log_x = False, log_y=False):
        """Plots subplots representing evolution of concentrations of desired compounds over time for different initial concentration of a reactant on the same figure to see 
        how initial concentration of a reactant affects the concentration profile of different compounds"""
        print(f"Plotting: {fig_name}...")
        fig, axes = self.initialize_fig_axes_objects(nrows, ncols, figsize)

        # Iterate through computed concentration dataframes and corresponding to them initial concentrations of studied reactant
        for i, (simulation_df, conc) in enumerate(zip(simulations_dfs, reactant_concentration_array)):
            ax = axes[i]

            # Get in shape of dataframe for plotting
            plot_df = pd.melt(simulation_df[compounds_for_plotting + ["time"]], id_vars = "time", var_name = "species", value_name = "conc")

            self.plot_axe_from_df_seaborn_groupplot(x = "time", 
                                                    y = "conc", 
                                                    xlabel = "time (h)",
                                                    ylabel = "concentration (M)", 
                                                    title = f"Initial C({reactant_to_study}) = {conc:.3e} M",
                                                    ax = ax,
                                                    palette = "muted",
                                                    hue = "species",
                                                    legend = True,
                                                    plot_df = plot_df,
                                                    log_y = log_y,
                                                    log_x = log_x
            )
        
        fig.suptitle(f"T = {temperature}K")

        self.hide_unused_subplots(fig, axes)
        return fig
    
    @staticmethod
    def calculate_R2_values(x_values, y_values, slope=None, intercept=None):
        """Calculates R2 value for given x and y values"""
        # Calculate slope if it wasn't calculated before
        if slope == None:
            slope, intercept = np.polyfit(x_values, y_values, 1)

        y_pred = PlotFunctions.get_points_for_trend_line(x_values, y_values, slope, intercept)

        ss_tot = np.sum((y_values - np.mean(y_values))**2)
        ss_res = np.sum((y_values - y_pred)**2)
        r2 = 1 - (ss_res / ss_tot)
        return r2
    
    @staticmethod 
    def get_points_for_trend_line(x_values, y_values, slope=None, intercept=None):
        """Calculates x and y points for a trend line for given x and y values"""
        x_values = np.asarray(x_values, dtype=float)
        y_values = np.asarray(y_values, dtype=float)
        
        if slope == None:
            slope, intercept = slope, intercept = np.polyfit(x_values, y_values, 1)

        trend_line = slope * x_values + intercept
        return trend_line

    def plot_ln_ri_vs_1_over_T(self, nrows, ncols, df, T_values_array, figsize, reactant_initial_concentration, reactant_to_study, fig_name):
        """Plots ln(ri) in a function of 1/T and a trendline for fluxes which R2 is greater than 0.9 (you can adjust it)"""
        print(f"Plotting: {fig_name}...")
        
        fig, axes = self.initialize_fig_axes_objects(nrows, ncols, figsize)

        df_at_reactant_c0 = df.loc[
                                    (df[f"c0({reactant_to_study})"] == reactant_initial_concentration) & 
                                    (df["R2"] > 0.9) & 
                                    (df["name"].str.contains(".Flux"))
                                  ].copy().reset_index()
        
        # print(df_at_reactant_c0.columns)

        for i, row in df_at_reactant_c0.iterrows():
            r2 = row["R2"]
            reaction = row["reaction"]
            slope = row["slope"]
            intercept = row["intercept"]
            Ea = row["Ea"] 

            # Taking and combining all 1/T columns to 1 np array 
            x_values = row[ [f"1/T({i})" for i in range(len(T_values_array))] ].to_numpy()

            # Taking and combining all ln values columns to 1 np array 
            y_values = row[ [f"ln_value({i})" for i in range(len(T_values_array))] ].to_numpy()

            intercept = row["intercept"]

            axes[i].scatter(x = x_values, 
                            y = y_values,
                            color = "red",
                            label = f"R2 = {r2:.4f}\nEa = {Ea:.2f}",
            )
            axes[i].set_title(reaction)
            axes[i].set_xlabel('1/T')
            axes[i].set_ylabel('ln(ri)')
            axes[i].legend()

            # print(x_values, y_values, slope, intercept)

            trend_line = self.get_points_for_trend_line(x_values, y_values, slope, intercept)

            axes[i].plot(x_values,
                        trend_line,
                        color = "blue"
            )

        fig.suptitle(f"c0({reactant_to_study}) = {reactant_initial_concentration}")
        self.hide_unused_subplots(fig, axes)
        return fig

    def plot_Ea_vs_reactant_c0(self, nrows, ncols, df, figsize, reactant_to_study, reactions_to_plot, fig_name, value_type, log_x=False):
        """Plots Ea in a function of c0 for specified fluxes/compound rates to see how initial concentration of a reactant affects the Ea of different steps in a system"""
        print(f"Plotting: {fig_name}...")
        fig, axes = self.initialize_fig_axes_objects(nrows, ncols, figsize)

        # Determine the grouping column and filter the dataframe
        group_col = "reaction" if value_type == "ri" else "compound"
        filtered_df = df[df[group_col].isin(reactions_to_plot)].groupby("name")

        for i, (name, group_df) in enumerate(filtered_df):
            
            # Extract plot data
            x_values = group_df[f"c0({reactant_to_study})"]
            y_values = group_df["Ea"]
            title = group_df[group_col].iloc[0]

            # Plot the data
            axes[i].plot(x_values, y_values, color = "red", marker = "o")
            axes[i].set_title(title)
            axes[i].set_xlabel(f"c0({reactant_to_study})")
            axes[i].set_ylabel("Ea")

            if log_x:
                axes[i].set_xscale("log")

        self.hide_unused_subplots(fig, axes)
        return fig
    
    def plot_drc_vs_c0(self, nrows, ncols, df, figsize, reactions, reactant_to_study, fig_name, temperature_value, log_x=False):
        print(f"Plotting: {fig_name}...")
        fig, axes = self.initialize_fig_axes_objects(nrows, ncols, figsize)

        filtered_df = df[df['temperature'] == temperature_value]
        c0_values = filtered_df[f"c0({reactant_to_study})"].values

        for i, reaction in enumerate(reactions):
            drc_values = filtered_df[reaction].values

            axes[i].plot(c0_values, drc_values, marker='o')

            axes[i].set_xlabel(f"c0({reactant_to_study})")
            axes[i].set_ylabel("drc_coefficients")
            axes[i].set_title(f'{reaction}')

            if log_x:
                axes[i].set_xscale("log")
        
        self.hide_unused_subplots(fig, axes)
        return fig

    @staticmethod
    def save_figure(fig_name):
        """Saves a figure in the dedicates directory"""
        try:
            print(f"Saving the figure: {fig_name}")
            plt.savefig(fig_name, dpi = 300)
            FileOperations.move_to_output_directory(IMAGES_SIMULATIONS_OUTPUT_DIR_NAME, fig_name)

        except Exception as e:
            print(f"Error saving figure {fig_name}: {e}")