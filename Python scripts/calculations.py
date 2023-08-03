
import pandas as pd
import numpy as np
from scipy import stats

def fetch_data(pSet):
    """
    Fetch the published and recomputed AUC values from a given PharmacoSet object.

    Parameters
    ----------
    pSet : object
        The pSet object from which to fetch data.

    Returns
    -------
    tuple of pd.DataFrame
        Two dataframes, the first containing published AUC values and the second containing recomputed AUC values.
    """

    df_pub = pSet.get_pub_AUC_values()
    df_recomp = pSet.recomp_Aucs
    return df_pub, df_recomp

class Correlation:

    

    @staticmethod
    def calculate_correlation(df1, df2, drug, pset1_name, pset2_name, method, type_):
        """
        Calculate the correlation between two datasets for a given drug using a specified method.

        Parameters
        ----------
        df1 : pd.DataFrame
            The first dataframe.    # Structered like pSet.recomp_Aucs 
        df2 : pd.DataFrame
            The second dataframe.
        drug : str
            The drug for which to calculate the correlation.
        pset1_name : str
            The name of the first pSet.
        pset2_name : str
            The name of the second pSet.
        method : str
            The correlation method to use ('spearman' or 'pearson').
        dataset_type : str
            The type of the dataset ('Published' or 'Recomputed').

        Returns
        -------
        tuple
            The calculated correlation coefficient and the sample size (n_cell).
        """
        df1_sub = df1[df1['drug'] == drug].drop(columns=['drug']).transpose()
        df1_sub.columns = [f'{type_}_{pset1_name}']

        df2_sub = df2[df2['drug'] == drug].drop(columns=['drug']).transpose()
        df2_sub.columns = [f'{type_}_{pset2_name}']

        merged_df = pd.concat([df1_sub, df2_sub], axis=1).dropna()
        
        size = len(merged_df)

        if method == "spearman":    
            correlation_coef = stats.spearmanr(merged_df[f'{type_}_{pset1_name}'], merged_df[f'{type_}_{pset2_name}'])
        elif method == "pearson":
            correlation_coef = stats.pearsonr(merged_df[f'{type_}_{pset1_name}'], merged_df[f'{type_}_{pset2_name}'])
        
        return correlation_coef, size

    @staticmethod
    def calculate_p_value(size_pub, size_recomp, correlation_coef_pub, correlation_coef_recomp):
        """
        Calculate the p-value for the significance of the difference between correlation coefficients of published and recomputed AUCs.

        Parameters
        ----------
        size_pub : int
            The sample size of the published AUC dataset (n_cell).
        size_recomp : int
            The sample size of the recomputed AUC dataset (n_cell).
        correlation_coef_pub : float
            The correlation coefficient for the published AUC.
        correlation_coef_recomp : float
            The correlation coefficient for the recomputed AUC.

        Returns
        -------
        float
            The calculated p-value.
        """
        # Fisher's z-transformations
        z1_pub = np.arctanh(correlation_coef_pub[0])
        z2_recomp = np.arctanh(correlation_coef_recomp[0])

        # Compute SE
        SE_diff = np.sqrt((1 / (size_pub - 3)) + (1 / (size_recomp - 3)))

        # Calculate Z-scores for the difference
        Z_score = (z1_pub - z2_recomp) / SE_diff

        p_value = (1 - stats.norm.cdf(abs(Z_score))) * 2
        
        return p_value
    

    @staticmethod
    def recomp_vs_pub_corr_table(pSet1, pSet2, method):
        """
        Create correlation table for all drugs between published and recomputed AUC values for two pSets.

        Parameters
        ----------
        pSet1 : object
            The first pSet object.
        pSet2 : object
            The second pSet object.
        method : str
            The correlation method to use ('spearman' or 'pearson').

        Returns
        -------
        pd.DataFrame
            DataFrame with correlation coefficients and p-values for each drug.
        """

        df_pub1, df_recomp1 = fetch_data(pSet1)
        df_pub2, df_recomp2 = fetch_data(pSet2)
        
        # Prepare results dataframe
        results = pd.DataFrame(columns=['drug', 'pub_corr', 'recomp_corr', 'pub_size', 'recomp_size', 'p-value'])
        
        for drug in df_recomp2['drug'].unique():
            try:
                correlation_coef_pub, size_pub = Correlation.calculate_correlation(df_pub1, df_pub2, drug, pSet1.name, pSet2.name, method, 'Published')
                correlation_coef_recomp, size_recomp = Correlation.calculate_correlation(df_recomp1, df_recomp2, drug, pSet1.name, pSet2.name, method, 'Recomputed')
                p_value = Correlation.calculate_p_value(size_pub, size_recomp, correlation_coef_pub, correlation_coef_recomp)
                
                new_row = pd.DataFrame({'drug': [drug], 
                            'pub_corr': [correlation_coef_pub[0]], 
                            'recomp_corr': [correlation_coef_recomp[0]], 
                            'pub_size': [size_pub], 
                            'recomp_size': [size_recomp], 
                            'p-value': [p_value]})

                # Append to results
                results = pd.concat([results, new_row], ignore_index=True)
            except:
                "Something wrong"
        
        return results

