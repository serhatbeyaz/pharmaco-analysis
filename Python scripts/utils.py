
import pandas as pd
import numpy as np
from scipy import integrate, stats

class Tools:
    
    @staticmethod
    def find_common_pairs(*PharmacoSets):
        """
        Finds the common cell line-drug pairs across multiple Pharmacological data sets (PharmacoSets).

        This method takes two or more PharmacoSets as input, and identifies cell line-drug pairs 
        that are present in all of them. It returns these common pairs as a set.

        Parameters
        ----------
        *PharmacoSets : 
            Variable length argument list of PharmacoSet objects. 

        Returns
        -------
        set
            A set containing tuples, where each tuple is a cell line-drug pair that is present 
            in all the input PharmacoSets.
        """
        # Initialize a set with cell_line-drug pairs from the first PharmacoSet
        common_pairs = set(PharmacoSets[0].data[['cell_line', 'drug']].itertuples(index=False))

        # Intersect with cell_line-drug pairs from each subsequent PharmacoSet
        for pset in PharmacoSets[1:]:
            pairs = set(pset.data[['cell_line', 'drug']].itertuples(index=False))
            common_pairs &= pairs

        return common_pairs

    @staticmethod
    def common_concentration_ranges(dataset1, dataset2) -> pd.DataFrame:
        """
        Calculate the common concentration ranges for drug-cell line pairs across two datasets.

        This method takes two datasets, combines them, and calculates the common minimum and 
        maximum concentrations (dose) for each drug-cell line pair. It returns a DataFrame containing 
        these common concentration ranges.

        Parameters
        ----------
        dataset1 : pd.DataFrame
            The first dataset to compare. It should have columns 'cell_line', 'drug', 'dose', and 'dataset'. 
            The 'dose' should be a positive number.

        dataset2 : pd.DataFrame
            The second dataset to compare. It should have columns 'cell_line', 'drug', 'dose', and 'dataset'. 
            The 'dose' should be a positive number.

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns 'cell_line_', 'drug_', 'common_min_dose', and 'common_max_dose', 
            where each row represents a unique combination of cell line and drug from the input datasets, and 
            'common_min_dose' and 'common_max_dose' are the common concentration range for that combination. 
            Rows with missing 'common_min_dose' or 'common_max_dose' values are excluded.
        """

        # Combine the two datasets
        combined_data = pd.concat([dataset1, dataset2])
        combined_data[['cell_line', 'drug']] = combined_data[['cell_line', 'drug']].apply(lambda x: x.str.lower())

        # Filter out rows with non-positive dose values
        combined_data = combined_data[combined_data['dose'] > 0]

        # Group data by 'cell_line', 'drug', and 'dataset' and calculate min and max dose for each group
        summary = combined_data.groupby(['cell_line', 'drug', 'dataset']).agg(min_dose=('dose', 'min'), 
                                                                               max_dose=('dose', 'max')).reset_index()

        # Pivot the dataset to have columns for min_dose and max_dose from both datasets
        pivot_table = pd.pivot_table(summary, values=['min_dose', 'max_dose'], index=['cell_line', 'drug'], 
                                     columns=['dataset']).reset_index()

        # Flatten the MultiIndex
        pivot_table.columns = ['_'.join(col).strip() for col in pivot_table.columns.values]

        dataset1_name = dataset1["dataset"].unique().item()
        dataset2_name = dataset2["dataset"].unique().item()

        # Calculate the common minimum and maximum dose
        pivot_table['common_min_dose'] = pivot_table[['min_dose_' + dataset1_name, 'min_dose_' + dataset2_name]].max(axis=1)
        pivot_table['common_max_dose'] = pivot_table[['max_dose_' + dataset1_name, 'max_dose_' + dataset2_name]].min(axis=1)

        # Filter out rows with missing common_min_dose or common_max_dose values
        pivot_table.dropna(subset=['common_min_dose', 'common_max_dose'], inplace=True)

        # Select the relevant columns
        common_dose_ranges = pivot_table[['cell_line_', 'drug_', 'common_min_dose', 'common_max_dose']]
        
        return common_dose_ranges
    
    @staticmethod
    def logistic_4PL(x, log_ec50, slope, front, back):
        model = (front - back) / (1 + 10 ** (slope * (x - log_ec50))) + back
        return model

    @staticmethod
    def compute_AUC(min_conc, max_conc, ec50, slope, front, back, conc_unit="uM") -> float:
        """
        Compute the Area Under the Curve (AUC) for a 4PL logistic regression model.

        This method takes the minimum and maximum concentrations (doses), the EC50, slope, front, and back 
        parameters of a 4PL logistic regression model, and calculates the AUC using the trapezoidal rule.
        The AUC is then normalized to a range of 0 to 1, and returned.

        Parameters
        ----------
        min_conc : float
            The minimum concentration value.

        max_conc : float
            The maximum concentration value.

        ec50 : float
            The EC50 value of the 4PL model.

        slope : float
            The slope of the 4PL model.

        front : float
            The front parameter of the 4PL model.

        back : float
            The back parameter of the 4PL model.

        conc_unit : str, optional
            The unit of the concentration values, by default "uM". 
            If set to "uM", the concentration values will be converted to M (by multiplying by 1e-6).

        Returns
        -------
        Optional[float]
            The normalized AUC of the 4PL model. If the minimum concentration is NaN, returns None.
        """
        # Convert the unit
        if conc_unit == "uM":
            min_conc = min_conc * 1e-06
            max_conc = max_conc * 1e-06

        if pd.isna(min_conc):
            return None

        # Concentration values
        min_concentration = np.log10(min_conc) 
        max_concentration = np.log10(max_conc)
        num_points = 1000  # Number of points to evaluate the function at

        # Generate a sequence of concentration points
        concentration_points = np.linspace(min_concentration, max_concentration, num_points)

        # Evaluate the 4PL model at the concentration points
        response_values = Tools.logistic_4PL(concentration_points, ec50, slope, front, back)

        # Calculate AUC using the trapezoidal rule
        auc = integrate.trapz(response_values, concentration_points)

        model_front = response_values[0]

        # Normalize AUC to the range of 0 to 1
        max_auc = model_front * (max_concentration - min_concentration)
        normalized_auc = auc / max_auc

        return 1 - normalized_auc
    






class DrugAnalysis:
    @staticmethod
    def compute_auc_across_cells(drug, PharmacoSet, common_dose_ranges, p_filter, log_fc) -> dict:
        """
        Computes the Area Under the Curve (AUC) for a given drug across multiple cell lines.

        This static method filters curve data, returned from the CurveCurator for the provided drug,
        based on a p-value and log_fc threshold and then calculates the AUC for each cell line within the 
        specified dose ranges taken from a DataFrame of common dose ranges.

        Parameters
        ----------
        drug : str
            The name of the drug for which AUC is to be calculated.
        PharmacoSet : PharmacoSet
            An instance of the PharmacoSet class containing pharmacological data.
        common_dose_ranges : pandas.DataFrame
            A DataFrame containing common dose ranges for each cell line-drug pair. This DataFrame 
            must include the following columns: 'cell_line_', 'drug_', 'common_min_dose', 'common_max_dose'.
        p_filter : float
            The p-value threshold for filtering curve data.
        log_fc : float
            The threshold for the logarithm of the fold change.

        Returns
        -------
        dict
            A dictionary where keys are cell line names and values are calculated AUCs. 
            If no curve data passing the p-value threshold is found, None is returned.

        Notes
        -----
        If the provided drug is not found in the common dose ranges for a given cell line, that cell line is skipped.
        """

        ## Get the Curve data
        curves = PharmacoSet.get_curve_data(drug).copy()

        
        if curves.empty:
            return None

        ## Convert p-value into Log p_value
        log_p = -np.log10(p_filter)

        curves = curves[(curves["Curve Log P_Value Corrected"] > log_p ) & (curves["Curve Fold Change"] < log_fc)]
        curves["Name"] = curves["Name"].str.lower()
  

        if curves.empty:
            print(f"No curve p-value < {p_filter} is found!")
            return None

        auc_list = {}

        for i in range(len(curves)):
            if drug in common_dose_ranges[common_dose_ranges['cell_line_'] == curves.iloc[i]["Name"]]['drug_'].values:

                conc_min = common_dose_ranges[(common_dose_ranges['cell_line_'] == curves.iloc[i]["Name"]) & 
                                              (common_dose_ranges['drug_'] == drug)]['common_min_dose'].values[0]

                conc_max = common_dose_ranges[(common_dose_ranges['cell_line_'] == curves.iloc[i]["Name"]) & 
                                              (common_dose_ranges['drug_'] == drug)]['common_max_dose'].values[0]

                params = curves.iloc[i][["Log EC50", "Curve Slope", "Curve Front", "Curve Back"]]

                auc = Tools.compute_AUC(min_conc=conc_min, max_conc=conc_max, ec50=params["Log EC50"],
                                        slope=params["Curve Slope"], front=params["Curve Front"],
                                        back=params["Curve Back"])

                auc_list[curves.iloc[i]["Name"]] = auc

            else:
                continue

        return auc_list
    
    @staticmethod
    def compute_median_absolute_deviation(auc_values):
        """
        Computes the median absolute deviation (MAD) of the given AUC values.

        This static method calculates and returns the MAD of a given dictionary of AUC values.

        Parameters
        ----------
        auc_values : dict
            A dictionary where keys are cell line names and values are calculated AUCs.

        Returns
        -------
        float
            The median absolute deviation of the AUC values.
        """
        # Extract the AUC values from the dictionary
        values = list(auc_values.values())

        # Compute and return the median absolute deviation
        return stats.median_abs_deviation(values)


