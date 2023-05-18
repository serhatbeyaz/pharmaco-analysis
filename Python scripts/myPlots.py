import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from pharmacoset import PharmacoSet
import pandas as pd
from scipy import stats
import seaborn as sns


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

        plt.title(f'Correlation of {drug} response in cell lines, {df1_name} vs. {df2_name} ')
        plt.text(0.1, 0.9, f'Spearman correlation: {spearsman_correlation[0]:.2f}', transform=ax.transAxes)
        plt.text(0.1, 0.85, f'Pearson correlation: {pearson_correlation[0]:.2f}', transform=ax.transAxes)

        plt.show()
        
    

        

    def barplot_corr_table(results):
    
        results_melt = results.melt(id_vars=['drug', 'p-value'], value_vars=['pub_corr', 'recomp_corr'], var_name='type', value_name='correlation')

        # Define a color palette with as many distinct colors as there are drugs
        n = len(results['drug'].unique())
        colors = sns.color_palette('hsv', n)

        fig, ax = plt.subplots(figsize=(8, 5))
        
        width = 0.35  # Width of the bars
        ind = np.arange(n)  # Location of the bars

        for i, drug in enumerate(results['drug'].unique()):
            data_pub = results_melt[(results_melt['drug'] == drug) & (results_melt['type'] == 'pub_corr')]
            data_recomp = results_melt[(results_melt['drug'] == drug) & (results_melt['type'] == 'recomp_corr')]
            
            ax.bar(ind[i] - width/2, data_pub['correlation'], width, color=colors[i], label='Published' if i == 0 else "")
            ax.bar(ind[i] + width/2, data_recomp['correlation'], width, hatch='///', fill=True, color=colors[i], label='Recomputed' if i == 0 else "", alpha=0.5, edgecolor='black')

        
        # Annotate bars with p-value
        for i, p in enumerate(results['p-value']):
            if p < 0.05:
                plt.text(i, results.loc[i, 'pub_corr'] if results.loc[i, 'pub_corr'] > results.loc[i, 'recomp_corr'] else results.loc[i, 'recomp_corr'], '*', ha='center')


        handles, labels = ax.get_legend_handles_labels()
        ax.set_title('Correlation Coefficients by Drug')
        ax.set_ylabel('Correlation Coefficient')
        ax.set_xlabel('Drug')
        ax.set_ylim(0, 1)  # Scale y-axis between 0 and 1
        ax.set_xticks(ind)
        ax.set_xticklabels(results['drug'].unique(), rotation=90)
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        

        fig.tight_layout()
        plt.show()




