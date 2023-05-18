import pandas as pd
import os
import numpy as np
from utils import DrugAnalysis


class PharmacoSet:
    """
    Class to manage pharmacological data sets.

    Attributes
    ----------
    name : str
        The name of the dataset, e.g. 'ccle', 'gdsc'.
    filepath : str
        The path to the dataset file. Should be processed data (Same input data used for CurveCurator)
    data : pd.DataFrame
        The data set as a pandas DataFrame.
    overlapping_data : pd.DataFrame or None
        Filtered data based on overlapping cell line-drug pairs.
        Initially None, filter_common_pairs() method can be used to assign.
    recomp_Aucs: pd.DataFrame or None
        Recomputed AUC values for all drugs across cell lines.
        Initially None, get_recomp_auc_values() method can be used to assign.
    """
    def __init__(self, name: str, filepath: str):
        """
        Initializes PharmacoSet with dataset name, file path, and loads the data into a DataFrame.

        Parameters
        ----------
        name : str
            The name of the dataset.
        filepath : str
            The path to the dataset file.
        """
        self.name = str(name)
        self.filepath = filepath
        self.data = pd.read_csv(filepath, sep= "\t")  # assuming a TSV file
        self.data.dropna(subset=['cell_line'], inplace=True)
        self.data['dataset'] = name
        self.overlapping_data = None
        self.recomp_Aucs = None
        
    def get_curve_data(self, drug: str) -> pd.DataFrame:
        """
        Loads the curve data outputed by CurveCurator for a specific drug.

        Parameters
        ----------
        drug : str
            The name of the drug.

        Returns
        -------
        pd.DataFrame
            Curve data for the specified drug.
            If the file is not found, returns an empty DataFrame.
        """
        path_to_curve = os.path.join("..", f"{self.name}_mapped_nCC", "curves", f"{drug}_curves.txt")
        try:
            curve_data = pd.read_csv(path_to_curve, sep= "\t")
            return curve_data
        except:
            FileNotFoundError 
            print(f"File is not found: {path_to_curve}")
            return pd.DataFrame()
    
    def get_pub_AUC_values(self) -> pd.DataFrame:
        """
        Fetches published AUC values.

        Returns
        -------
        pd.DataFrame
            Published AUC values for the PharmacoSet.
        """
        path_to_file = os.path.join("..", "Auc_values_published", f"AUC_values_{self.name}_published.csv")
        auc_values = pd.read_csv(path_to_file)
        return auc_values
        
        
    def filter_common_pairs(self, common_pairs: set) -> None:
        """
        Filters the data to include only cell line-drug pairs that are common across multiple data sets.

        Parameters
        ----------
        common_pairs : set
            Set of common cell line-drug pairs.
            Tools.find_common_pairs() method can be used to create required set.

        Returns
        -------
        None
            This method doesn't return anything, but overwrites the `overlapping_data` attribute of the class instance.
        """
        common_pairs_set = set(common_pairs)
        # Create a mask where each element is True if the corresponding cell line-drug pair is in common_pairs
        mask = self.data.apply(lambda row: (row['cell_line'], row['drug']) in common_pairs_set, axis=1)
        # Filter the DataFrame using the mask
        self.overlapping_data = self.data[mask]
    
    def get_recomp_auc_values(self, common_concentration, p_filter = 0.05, log_fc = -0.05):
        """
        Computes AUC values for a list of drugs using the parameters from CurveCurator.

        Parameters
        ----------
        common_concentration : pd.DataFrame
            DataFrame of common concentrations.
            Can be generated by using Tools.common_concentration_ranges().
        p_filter : float, optional
            p_value threshold to filter the curves
        log_fc : float, optional
            Log Fold Change threshold to filter out the curves

        Returns
        -------
        pd.DataFrame
            Re-computed AUC values for the specified drugs.
        """
        drug_list = self.overlapping_data['drug'].unique()
        all_drugs = {}

        for drug in drug_list:
            auc_drug = DrugAnalysis.compute_auc_across_cells(drug= drug, PharmacoSet=self, common_dose_ranges=common_concentration, p_filter= p_filter, log_fc= log_fc)
            if auc_drug == None:
                continue
            all_drugs[drug] = auc_drug

        df = pd.DataFrame(all_drugs).T
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'drug'}, inplace=True)

        self.recomp_Aucs = df
        return df
        
        


