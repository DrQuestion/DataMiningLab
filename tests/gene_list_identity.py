import os
import pandas as pd
root = "C:/Users/aless/Documents/DataMiningLab/gdc_download_20191010_140437.280422"
prec_series = pd.Series()
sample_no = 0
for folder in os.listdir(root):
    files = os.listdir(os.path.join(root, folder))
    for file in files:
        file_path = os.path.join(root, folder, file)
        if file.endswith('.FPKM.txt'):
            sample_expr = pd.read_csv(file_path, sep='\t')
            # to look if all genes_ID are the sames
            if sample_no < 1:
                prec_series = sample_expr.iloc[:, 0]
                sample_no = 1
            else:
                series = sample_expr.iloc[:, 0]
                if not prec_series.equals(series):
                    print(file)
                prec_series = series
