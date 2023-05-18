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
        