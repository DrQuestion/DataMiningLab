import os
import pandas as pd
expression_data_path = r"C:/Users/aless/Documents/raw_counts"
sample_sheet_path = r"C:/Users/aless/Documents/gdc_sample_sheet.2019-10-23.tsv"
matrix = pd.DataFrame()
sample_no = 0
file_to_barcode = dict()
patients = set()
patients_with_replicas = set()
sheet_table = pd.read_csv(sample_sheet_path, sep="\t")
for row_id, row in sheet_table.iterrows():
    barcode = row["Sample ID"]
    file_to_barcode[row["File Name"][:-3]] = barcode
    patient = barcode[8:-4]
    if patient not in patients:
        patients.add(barcode[8:-4])
    else:
        patients_with_replicas.add(patient)
for folder in os.listdir(expression_data_path):
    files = os.listdir(os.path.join(expression_data_path, folder))
    for file in files:
        file_path = os.path.join(expression_data_path, folder, file)
        if file.endswith('.counts'):
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
print(matrix.shape)
matrix.to_csv(r"C:/Users/aless/Documents/DataMiningLab/rawCountsExpressionMatrix.tsv", sep='\t')
