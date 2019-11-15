import os
import numpy as np
import pandas as pd
expression_data_path = r"C:/Users/aless/Documents/DataMiningLab/gdc_download_20191010_140437.280422"
sample_sheet_path = r"C:/Users/aless/Documents/DataMiningLab/gdc_sample_sheet.2019-10-11.tsv"


def mat_gen(expression_matrix_path=expression_data_path, sample_sheet_path=sample_sheet_path, data_type='fpkms'):
    matrix = pd.DataFrame()
    sample_no = 0
    file_to_barcode = dict()
    patients = set()
    patients_with_replicas = set()
    sheet_table = pd.read_csv(sample_sheet_path, sep="\t")
    extension = '.FPKM.txt'
    if data_type == 'raw_counts':
        extension = '.counts'
    for row_id, row in sheet_table.iterrows():
        barcode = row["Sample ID"]
        file_to_barcode[row["File Name"][:-3]] = barcode
        patient = barcode[8:-4]
        if patient not in patients:
            patients.add(barcode[8:-4])
        else:
            patients_with_replicas.add(patient)
    for folder in os.listdir(expression_matrix_path):
        files = os.listdir(os.path.join(expression_matrix_path, folder))
        for file in files:
            file_path = os.path.join(expression_matrix_path, folder, file)
            if file.endswith(extension):
                barcode = file_to_barcode[file]
                if barcode[8:-4] not in patients_with_replicas:
                    sample_expr = pd.read_csv(file_path, sep='\t')
                    if sample_no < 1:
                        matrix['Gene_ID'] = sample_expr.iloc[:, 0]
                        matrix[barcode] = sample_expr.iloc[:, 1]
                        sample_no = 1
                    else:
                        matrix[barcode] = sample_expr.iloc[:, 1]
                else:
                    continue
    matrix = matrix.set_index('Gene_ID')
    matrix.to_csv(r"C:/Users/aless/Documents/DataMiningLab/expressionMatrix.tsv", sep='\t')
    return matrix


def matrix_trim(matrix, zero_ratio=1):
    if zero_ratio == 1:
        matrix = matrix.loc[~(matrix == 0).all(axis=1)]
        return matrix
    else:
        more_than_ratio = list()
        ncols = matrix.shape[0]
        for i, row in matrix.iterrows():
            ratio = np.count_nonzero(row == 0)/ncols
            if ratio >= zero_ratio:
                more_than_ratio.append(i)
        matrix = matrix.drop(index=more_than_ratio)
    return matrix
