import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from pharmacoset import PharmacoSet
import pandas as pd
from scipy import stats
import seaborn as sns
import calculations
import utils
import itertools


class Plotting:

    @staticmethod
    def correlation_plot(pSet1: PharmacoSet,  pSet2: PharmacoSet, drug:str, method = 'recomputed'):   # method should be either 'recomputed' or 'published'
        """
        Compare dataframes to evaluate drug response correlations between two pharmacological datasets.

        This method generates a scatterplot comparing drug responses in cell lines between two PharmacoSet objects.
        The method also computes and displays Spearman and Pearson correlation coefficients on the plot.

        Parameters
        ----------
        pSet1 : PharmacoSet
            First pharmacological dataset.
        pSet2 : PharmacoSet
            Second pharmacological dataset.
        drug : str
            Name of the drug to compare.
        method : str, optional
            Specifies the method for comparison. It can be either 'recomputed' (default) or 'published'.

        Returns
        -------
        None. Displays a scatterplot comparing the two datasets, with data points labeled according to their response (sensitive or resistant). AUC > 0.2 == sensitive 
        The plot also displays lines indicating the thresholds at x=0.2, y=0.2, and x=y, as well as Spearman and Pearson correlation coefficients.

        Notes
        -----
        The 'published' method uses the get_pub_AUC_values method of the PharmacoSet objects,
        while the 'recomputed' method uses the recomp_Aucs attribute.
        """

        df1_name = pSet1.name.upper()
        df2_name = pSet2.name.upper()

        if method == 'published':
            df1 = pSet1.get_pub_AUC_values()
            df2 = pSet2.get_pub_AUC_values()
        else:
            df1 = pSet1.recomp_Aucs
            df2 = pSet2.recomp_Aucs
            
        df1_drug = df1[df1['drug'] == drug].drop(columns=['drug']).transpose()
        df1_drug.columns = [f'Recomputed_{df1_name}']

        df2_drug = df2[df2['drug'] == drug].drop(columns=['drug']).transpose()
        df2_drug.columns = [f'Recomputed_{df2_name}']

        merged_df = pd.concat([df1_drug, df2_drug], axis=1).dropna()
        
        merged_df['Label'] = np.select(
            [
                (merged_df[f'Recomputed_{df1_name}'] <= 0.2) & (merged_df[f'Recomputed_{df2_name}'] <= 0.2),
                (merged_df[f'Recomputed_{df1_name}'] > 0.2) & (merged_df[f'Recomputed_{df2_name}'] > 0.2),
                (merged_df[f'Recomputed_{df1_name}'] > 0.2) & (merged_df[f'Recomputed_{df2_name}'] <= 0.2),
                (merged_df[f'Recomputed_{df1_name}'] <= 0.2) & (merged_df[f'Recomputed_{df2_name}'] > 0.2)
            ], 
            [
                "Both resistant",
                "Both sensitive",
                f"{df1_name} sensitive, {df2_name} resistant",
                f"{df1_name} resistant, {df2_name} sensitive"
            ], 
            default=np.nan
        )

        spearsman_correlation = stats.spearmanr(merged_df[f'Recomputed_{df1_name}'], merged_df[f'Recomputed_{df2_name}'])
        pearson_correlation = stats.pearsonr(merged_df[f'Recomputed_{df1_name}'], merged_df[f'Recomputed_{df2_name}'])

        labels_sorted = sorted(merged_df['Label'].unique())

        fig, ax = plt.subplots(figsize=(10,8))

        sns.scatterplot(data=merged_df, x=f'Recomputed_{df1_name}', y=f'Recomputed_{df2_name}', hue='Label', hue_order=labels_sorted, ax=ax)

        ax.axline((0, 0), slope=1, color='r', linestyle='--')  # x=y line
        ax.axhline(0.2, color='b', linestyle='--')  # horizontal line at y=0.2
        ax.axvline(0.2, color='g', linestyle='--')  # vertical line at x=0.2

        ax.set_xlim([0, 1]) # set the x axis limit
        ax.set_ylim([0, 1]) # set the y axis limit

        plt.xlabel(f'{method}_{df1_name.upper()}')
        plt.ylabel(f'{method}_{df2_name.upper()}')

        plt.title(f'Correlation of {drug} response in cell lines, {df1_name} vs. {df2_name} ')
        plt.text(0.1, 0.9, f'Spearman correlation: {spearsman_correlation[0]:.2f}', transform=ax.transAxes)
        plt.text(0.1, 0.85, f'Pearson correlation: {pearson_correlation[0]:.2f}', transform=ax.transAxes)

        plt.show()
        
    


    def barplot_corr_table(results, method):
        # Sort results by 'recomp_corr' in descending order and reset index
        results = results.sort_values(by='recomp_corr', ascending=False).reset_index(drop=True)

        results_melt = results.melt(id_vars=['drug', 'p-value'], value_vars=['pub_corr', 'recomp_corr'], var_name='type', value_name='correlation')

        fig, ax = plt.subplots(figsize=(16, 10))
        
        width = 0.35  # Width of the bars
        ind = np.arange(len(results['drug'].unique()))  # Location of the bars

        pub_color = "skyblue"
        recomp_color = "orange"

        for i in range(len(results['drug'].unique())):
            data_pub = results_melt[(results_melt['drug'] == results['drug'].iloc[i]) & (results_melt['type'] == 'pub_corr')]
            data_recomp = results_melt[(results_melt['drug'] == results['drug'].iloc[i]) & (results_melt['type'] == 'recomp_corr')]
            
            ax.bar(ind[i] - width/2, data_pub['correlation'], width, color=pub_color, label='Published' if i == 0 else "")
            ax.bar(ind[i] + width/2, data_recomp['correlation'], width, hatch='///', fill=True, color=recomp_color, label='Recomputed' if i == 0 else "", edgecolor='black')

        # # Annotate bars with p-value
        # for i, p in enumerate(results['p-value']):
        #     if p < 0.05:
        #         pub_height = results.loc[i, 'pub_corr']
        #         recomp_height = results.loc[i, 'recomp_corr']
        #         plt.text(i, max(pub_height, recomp_height), '*', ha='center')  # No offset, place asterisk at top of higher bar

        handles, labels = ax.get_legend_handles_labels()
        ax.set_title(f'{method} Correlation Coefficients by Drug')
        ax.set_ylabel('Correlation Coefficient')
        ax.set_xlabel('Drug')
        ax.set_ylim(0, 1)  # Scale y-axis between 0 and 1
        ax.set_xticks(ind)
        ax.set_xticklabels(results['drug'], rotation=90)
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        
        fig.tight_layout()
        plt.show()


    def plot_AUC_dist_pub_vs_recomp(pSet):
        # Fetch the data
        df_pub, df_recomp = calculations.fetch_data(pSet)

        # Convert cell line names to lowercase
        df_pub.columns = map(str.lower, df_pub.columns)
        df_recomp.columns = map(str.lower, df_recomp.columns)

        # Reshape both dataframes to long format
        df_pub_long = df_pub.melt(id_vars='drug', var_name='cell_line', value_name='AUC_pub')
        df_recomp_long = df_recomp.melt(id_vars='drug', var_name='cell_line', value_name='AUC_recomp')

        # Merge both datasets based on 'drug' and 'cell_line'
        df_merged = pd.merge(df_pub_long, df_recomp_long, on=['drug', 'cell_line'])

        # Drop rows with NA values in either 'AUC_pub' or 'AUC_recomp'
        df_merged.dropna(subset=['AUC_pub', 'AUC_recomp'], inplace=True)

        # Unique drugs
        drugs = df_merged['drug'].unique()
        drugs.sort()

        if len(drugs) == 0:
            print("No unique drugs found in the merged dataframe.")
            return

        # Find global min and max AUC values to set the x-axis range for histograms
        global_min = df_merged[['AUC_pub', 'AUC_recomp']].min().min()
        global_max = df_merged[['AUC_pub', 'AUC_recomp']].max().max()

        # Determine bin edges based on global min/max and desired number of bins
        bin_edges = np.linspace(global_min, global_max, 21)

        # Plot
        num_plots = len(drugs)
        num_cols = 4  # We will arrange the plots in 4 columns.
        num_rows = num_plots // num_cols if num_plots % num_cols == 0 else num_plots // num_cols + 1
        
        fig, axs = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))

        for i, drug in enumerate(drugs):
            ax = axs[i // num_cols, i % num_cols]
            sns.histplot(data=df_merged[df_merged['drug'] == drug], x='AUC_pub', kde=True, ax=ax, label='Published', color='blue', alpha=0.5, bins=bin_edges)
            sns.histplot(data=df_merged[df_merged['drug'] == drug], x='AUC_recomp', kde=True, ax=ax, label='Recomputed', color='red', alpha=0.5, bins=bin_edges)
            ax.set_title(drug)
            ax.legend()

        # If there are less plots than slots in the grid, remove the empty plots
        if num_plots % num_cols != 0:
            for idx in range(num_plots, num_rows * num_cols):
                fig.delaxes(axs.flatten()[idx])

        plt.tight_layout()
        plt.show()

    

    def plot_waterfall(df, drug_name):
        # Melt the dataframe to a long format
        df_melted = df.melt(id_vars='drug', var_name='cell_line', value_name='AUC')

        # Drop rows with NaN AUC values
        df_melted.dropna(subset=['AUC'], inplace=True)
        
        # Filter for the specified drug
        df_drug = df_melted[df_melted['drug'] == drug_name]
        
        # Sort values by AUC for the waterfall effect
        df_drug = df_drug.sort_values(by='AUC')
        
        # Create the plot
        plt.figure(figsize=(15, 6))  
        plt.bar(df_drug['cell_line'], df_drug['AUC'], color='skyblue')
        plt.xticks(rotation=90, fontsize=8) 
        plt.title(f"Waterfall plot for {drug_name}")
        plt.ylim(0,1)
        plt.ylabel("AUC")
        plt.tight_layout()  
        plt.show()



    def plot_sigmoid_curves(pSet1, pSet2, drug, cell_line, common_conc):

        curve1 = pSet1.get_curve_data(drug)
        curve2 = pSet2.get_curve_data(drug)

        curve1 = curve1[curve1["Name"] == cell_line]
        curve2 = curve2[curve2["Name"] == cell_line]


        curve1_doses = pSet1.overlapping_data[(pSet1.overlapping_data['drug'] == drug) & (pSet1.overlapping_data['cell_line'] == cell_line)]['dose']

        org_len = len(curve1_doses)

        curve1_doses = np.array(list(itertools.dropwhile(lambda x: x == 0, curve1_doses)))

        removed_elements = org_len - len(curve1_doses)

        curve1_response = pSet1.overlapping_data[(pSet1.overlapping_data['drug'] == drug) & (pSet1.overlapping_data['cell_line'] == cell_line)]['response'].iloc[removed_elements:]


        curve2_doses = pSet2.overlapping_data[(pSet2.overlapping_data['drug'] == drug) & (pSet2.overlapping_data['cell_line'] == cell_line)]['dose']

        org_len2 = len(curve2_doses)    

        curve2_doses = np.array(list(itertools.dropwhile(lambda x: x == 0, curve2_doses)))

        removed_elements2 = org_len2 - len(curve2_doses)

        curve2_response = pSet2.overlapping_data[(pSet2.overlapping_data['drug'] == drug) & (pSet2.overlapping_data['cell_line'] == cell_line)]['response'].iloc[removed_elements2:]
        

        max_conc1 = np.log10(curve1_doses.max() * 1e-6)
        concentration_points1 = np.linspace(np.log10(curve1_doses.min() * 2e-7), max_conc1, 400)

        # Fetch parameters from dataset1 and plot curve
        log_ec50_1, slope_1, front_1, back_1 = curve1['Log EC50'].item(), curve1['Curve Slope'].item(), curve1['Curve Front'].item(), curve1['Curve Back'].item()

        plt.figure(figsize=(7,5))

        y1 = utils.Tools.logistic_4PL(concentration_points1, log_ec50_1, slope_1, front_1, back_1)
        plt.plot(10 ** concentration_points1, y1, label= f'{pSet1.name}')
        plt.scatter(curve1_doses * 1e-6, curve1_response, s= 15)


        max_conc2 = np.log10(curve2_doses.max() * 1e-6)
        concentration_points2 = np.linspace(np.log10(curve2_doses.min() * 2e-7), max_conc2, 400)

        # Fetch parameters from dataset2 and plot curve
        log_ec50_2, slope_2, front_2, back_2 = curve2['Log EC50'].item(), curve2['Curve Slope'].item(), curve2['Curve Front'].item(), curve2['Curve Back'].item()
        y2 = utils.Tools.logistic_4PL(concentration_points2, log_ec50_2, slope_2, front_2, back_2)

        plt.plot(10 ** concentration_points2, y2, label=f'{pSet2.name}')
        plt.scatter(curve2_doses * 1e-6, curve2_response, s = 15)


        common_conc = common_conc[(common_conc["cell_line_"] == cell_line.lower()) & (common_conc["drug_"] == drug)]
        common_min, common_max = common_conc['common_min_dose'].item() * 1e-6, common_conc['common_max_dose'].item() * 1e-6

        
        plt.title(f'{drug} - {cell_line} Curves')
        plt.xlabel('Concentration')
        plt.ylabel('Response')
        plt.xscale('log')
        plt.ylim(0, 2)
        plt.axvline(x = common_min, color = 'b', label = 'common_min', linestyle='--')
        plt.axvline(x = common_max, color = 'r', label = 'common_max', linestyle='--')
        plt.legend()
        plt.grid(True)
        plt.show()


    def multiple_sigmoid_curves(pSet1, pSet2, drug, cell_line, common_conc, ax):
        curve1 = pSet1.get_curve_data(drug)
        curve2 = pSet2.get_curve_data(drug)

        curve1 = curve1[curve1["Name"] == cell_line]
        curve2 = curve2[curve2["Name"] == cell_line]


        curve1_doses = pSet1.overlapping_data[(pSet1.overlapping_data['drug'] == drug) & (pSet1.overlapping_data['cell_line'] == cell_line)]['dose'].iloc[1:]
        curve1_response = pSet1.overlapping_data[(pSet1.overlapping_data['drug'] == drug) & (pSet1.overlapping_data['cell_line'] == cell_line)]['response'].iloc[1:]

        curve2_doses = pSet2.overlapping_data[(pSet2.overlapping_data['drug'] == drug) & (pSet2.overlapping_data['cell_line'] == cell_line)]['dose'].iloc[1:]
        curve2_response = pSet2.overlapping_data[(pSet2.overlapping_data['drug'] == drug) & (pSet2.overlapping_data['cell_line'] == cell_line)]['response'].iloc[1:]

        max_conc1 = np.log10(curve1_doses.max() * 1e-6)
        concentration_points1 = np.linspace(np.log10(curve1_doses.min() * 2e-7), max_conc1, 400)

        # Fetch parameters from dataset1 and plot curve
        log_ec50_1, slope_1, front_1, back_1 = curve1['Log EC50'].item(), curve1['Curve Slope'].item(), curve1['Curve Front'].item(), curve1['Curve Back'].item()


        y1 = utils.Tools.logistic_4PL(concentration_points1, log_ec50_1, slope_1, front_1, back_1)
        ax.plot(10 ** concentration_points1, y1, label= f'{pSet1.name}')
        ax.scatter(curve1_doses * 1e-6, curve1_response, s= 15)

        max_conc2 = np.log10(curve2_doses.max() * 1e-6)
        concentration_points2 = np.linspace(np.log10(curve2_doses.min() * 2e-7), max_conc2, 400)

        # Fetch parameters from dataset2 and plot curve
        log_ec50_2, slope_2, front_2, back_2 = curve2['Log EC50'].item(), curve2['Curve Slope'].item(), curve2['Curve Front'].item(), curve2['Curve Back'].item()
        y2 = utils.Tools.logistic_4PL(concentration_points2, log_ec50_2, slope_2, front_2, back_2)

        ax.plot(10 ** concentration_points2, y2, label=f'{pSet2.name}')
        ax.scatter(curve2_doses * 1e-6, curve2_response, s = 15)

        common_conc = common_conc[(common_conc["cell_line_"] == cell_line.lower()) & (common_conc["drug_"] == drug)]
        common_min, common_max = common_conc['common_min_dose'].item() * 1e-6, common_conc['common_max_dose'].item() * 1e-6

        ax.set_title(f'{drug} - {cell_line} Curves')
        ax.set_xlabel('Concentration')
        ax.set_ylabel('Response')
        ax.set_xscale('log')
        ax.set_ylim(0, 2)
        ax.axvline(x = common_min, color = 'b', label = 'common_min', linestyle='--')
        ax.axvline(x = common_max, color = 'r', label = 'common_max', linestyle='--')
        ax.legend()
        ax.grid(True)


    
    @staticmethod
    def correlation_plot_whole_pSet(pSet1,  pSet2, method = 'recomputed', drugs_to_exclude = None):  
        """
        Compare dataframes to evaluate response correlations between two pharmacological datasets.

        This method generates a scatterplot comparing responses in cell lines between two PharmacoSet objects.
        The method also computes and displays Spearman and Pearson correlation coefficients on the plot.

        Parameters
        ----------
        pSet1 : PharmacoSet
            First pharmacological dataset.
        pSet2 : PharmacoSet
            Second pharmacological dataset.
        method : str, optional
            Specifies the method for comparison. It can be either 'recomputed' (default) or 'published'.
        drugs_to_exclude : list, optional
            List of drug names to exclude from comparison.

        Returns
        -------
        None. Displays a scatterplot comparing the two datasets as well as Spearman and Pearson correlation coefficients.
        """

        df1_name = pSet1.name.upper()
        df2_name = pSet2.name.upper()

        if method == 'published':
            df1 = pSet1.get_pub_AUC_values()
            df2 = pSet2.get_pub_AUC_values()
        else:
            df1 = pSet1.recomp_Aucs
            df2 = pSet2.recomp_Aucs
                
        df1 = df1.set_index('drug').stack().reset_index().rename(columns={'level_1': 'cell_line', 0: f'Recomputed_{df1_name}'})
        df2 = df2.set_index('drug').stack().reset_index().rename(columns={'level_1': 'cell_line', 0: f'Recomputed_{df2_name}'})

        if drugs_to_exclude is not None:
            df1 = df1[~df1['drug'].isin(drugs_to_exclude)]
            df2 = df2[~df2['drug'].isin(drugs_to_exclude)]

        merged_df = pd.merge(df1, df2, on=['drug', 'cell_line'])

        spearsman_correlation = stats.spearmanr(merged_df[f'Recomputed_{df1_name}'], merged_df[f'Recomputed_{df2_name}'])
        pearson_correlation = stats.pearsonr(merged_df[f'Recomputed_{df1_name}'], merged_df[f'Recomputed_{df2_name}'])

        fig, ax = plt.subplots(figsize=(10,8))

        sns.scatterplot(data=merged_df, x=f'Recomputed_{df1_name}', y=f'Recomputed_{df2_name}', ax=ax)

        ax.axline((0, 0), slope=1, color='r', linestyle='--')  # x=y line
        ax.axhline(0.2, color='b', linestyle='--')  # horizontal line at y=0.2
        ax.axvline(0.2, color='g', linestyle='--')  # vertical line at x=0.2

        ax.set_xlim([0, 1]) 
        ax.set_ylim([0, 1]) 

        plt.xlabel(f'{method}_{df1_name.upper()}')
        plt.ylabel(f'{method}_{df2_name.upper()}')

        plt.title(f'Overall Correlation of response in cell lines, {df1_name} vs. {df2_name} ')
        plt.text(0.1, 0.9, f'Spearman correlation: {spearsman_correlation[0]:.2f}', transform=ax.transAxes)
        plt.text(0.1, 0.85, f'Pearson correlation: {pearson_correlation[0]:.2f}', transform=ax.transAxes)

        plt.show()
