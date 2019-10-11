import os
import pandas as pd
root="C:/Users/aless/Documents/DataMiningLab/gdc_download_20191010_140437.280422"
sample_sheet_path = "C:/Users/aless/Documents/DataMiningLab/gdc_sample_sheet.2019-10-11"
matrix=pd.DataFrame()
sample_no=0
file_prec = ""
file_to_barcode = dict()

for folder in os.listdir(root):
    files=os.listdir(os.path.join(root, folder))
    for file in files:
        barcode = sheet
        file_path = os.path.join(root, folder, file)
        if file.endswith('.FPKM.txt'):
            sample_expr=pd.read_csv(file_path, sep='\t')
            if sample_no < 1:
                matrix['Gene_ID']= sample_expr.iloc[:,0]
                matrix[file] = sample_expr.iloc[:,1]
                sample_no = 1
            else:
                matrix[file] = sample_expr.iloc[:,1]
