import os
import pandas as pd
root=r"C:/Users/aless/Documents/DataMiningLab/gdc_download_20191010_140437.280422"
sample_sheet_path = r"C:/Users/aless/Documents/DataMiningLab/gdc_sample_sheet.2019-10-11.tsv"
matrix=pd.DataFrame()
sample_no=0
file_prec = ""
file_to_barcode = dict()
sheet_table = pd.read_csv(sample_sheet_path, sep="\t")
print(sheet_table.shape)
for id, row in sheet_table.iterrows():
    file_to_barcode[row["File Name"][:-3]] = row["Sample ID"]
for file1 in file_to_barcode:
    for file2 in file_to_barcode:
        if file1 != file2 and file_to_barcode[file1] == file_to_barcode[file2]:
            print(file1, file2)
cnt=0
for folder in os.listdir(root):
    files=os.listdir(os.path.join(root, folder))
    flag = 0
    for file in files:
        file_path = os.path.join(root, folder, file)

        if file.endswith('.FPKM.txt'):
            flag = 1
            cnt+=1
            barcode = file_to_barcode[file]
            sample_expr=pd.read_csv(file_path, sep='\t')
            if sample_no < 1:
                matrix['Gene_ID']= sample_expr.iloc[:,0]
                matrix[barcode] = sample_expr.iloc[:,1]
                sample_no = 1
            else:
                matrix[barcode] = sample_expr.iloc[:,1]
    if flag == 0:
        print(folder)
print(matrix)