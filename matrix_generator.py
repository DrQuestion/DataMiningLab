import os
import numpy as np
import pandas as pd
expression_data_path = r"C:/Users/aless/Documents/UniTn/DataMiningLab/gdc_download_20191010_140437.280422"
sample_sheet_path = r"C:/Users/aless/Documents/UniTn/DataMiningLab/gdc_sample_sheet.2019-10-11.tsv"
matrix_path = r"C:/Users/aless/Documents/UniTn/DataMiningLab/expressionMatrix.tsv"
FINAL_NUMBER_OF_FEATURES = 7000
INDIVUMED_TSS = ['AA', 'AG']
NO_METADATA = ['TCGA-5M-AATA', 'TCGA-5M-AAT5', 'TCGA-F5-6810']


def matrix_generator(expression_matrix_path=expression_data_path, sample_sheet_path=sample_sheet_path,
                     barcodes_filter=[], data_type='fpkms', remove_indivumed=True, save_matrix=True):
    if not remove_indivumed:
        INDIVUMED_TSS = []
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

                if barcodes_filter:
                    if barcode[8:-4] not in patients_with_replicas and barcode[:-4] in barcodes_filter:
                        sample_expr = pd.read_csv(file_path, sep='\t')
                        if sample_no < 1:
                            matrix['Gene_ID'] = sample_expr.iloc[:, 0]
                            matrix[barcode] = sample_expr.iloc[:, 1]
                            sample_no = 1
                        else:
                            matrix[barcode] = sample_expr.iloc[:, 1]

                else:
                    if barcode[8:-4] not in patients_with_replicas and barcode[5:7] not in INDIVUMED_TSS \
                            and barcode[:-4] not in NO_METADATA:
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
    if save_matrix:
        matrix.to_csv(r"C:/Users/aless/Documents/UniTn/DataMiningLab/expressionMatrix.tsv", sep='\t')
    return matrix


def filter_generator(path_objective):
    objective = pd.read_csv(path_objective, sep='\t', header=0)
    barcodes = list(objective.iloc[:, 1])
    return barcodes


def matrix_trim(matrix, zero_ratio=0.9):
    if zero_ratio == 1:
        matrix = matrix.loc[~(matrix == 0).all(axis=1)]
        return matrix
    else:
        more_than_ratio = list()
        ncols = matrix.shape[1]
        for i, row in matrix.iterrows():
            count = np.count_nonzero(row == 0)
            ratio = count/ncols
            if ratio >= zero_ratio:
                more_than_ratio.append(i)
        matrix = matrix.drop(index=more_than_ratio)

    print(f'Trimming matrix for zero_ratio={zero_ratio} results in a {matrix.shape} shaped matrix')
    return matrix


def matrix_normalization(matrix, norm_method='log'):
    if norm_method == 'deseq2':
        matrix = deseq2_norm(matrix)
    elif norm_method == 'log':
        matrix = log_norm(matrix)
    elif norm_method == 'tpm':
        matrix = fpkm_to_tpm_norm(matrix)
    elif norm_method == 'no_norm':
        pass
    else:
        print('Wrong Norm')
        exit()
    return matrix


def log_norm(matrix):
    # log normalizes a table
    return np.log(matrix + 0.00000001)


def fpkm_to_tpm_norm(matrix):
    # converts a (samples x genes) matrix of fpkms to tpms
    total_fpkms = list()
    for i, row in matrix.iterrows():
        tot_sample_fpkm = sum(row)
        total_fpkms.append(tot_sample_fpkm)
    matrix = matrix.div(total_fpkms, axis='rows')
    return matrix * (10 ** 6)


def deseq2_norm(matrix):
    # converts a row counts matrix to a DESeq2 normalized matrix
    mirror_matrix = np.log(np.float16(matrix.T))
    mirror_matrix = pd.DataFrame(data=mirror_matrix, index=matrix.T.index.values, columns=matrix.T.columns.values)
    any_zero_index = list()
    for i, row in mirror_matrix.iterrows():
        if row.isin([-np.inf]).any():
            any_zero_index.append(i)
        else:
            geo_mean = sum(row)/len(row)
            mirror_matrix.loc[i, :] = row-geo_mean
    mirror_matrix = mirror_matrix.drop(index=any_zero_index)
    scaling_factors = np.exp(mirror_matrix.median())
    return matrix.div(scaling_factors, axis='rows')


def matrix_load(path=matrix_path):
    matrix = pd.read_csv(path, sep='\t', index_col=0, header=0)
    return matrix


def filter_by_variance(matrix, final_number_of_features=FINAL_NUMBER_OF_FEATURES):
    variances = list()
    minor_variance = list()
    for i, row in matrix.iterrows():
        variances.append(row.var())
    threshold_var = sorted(variances)[-(final_number_of_features+1)]
    for i, row in matrix.iterrows():
        if row.var() <= threshold_var:
            minor_variance.append(i)
    matrix = matrix.drop(index=minor_variance)
    matrix.to_csv(rf"C:/Users/aless/Documents/UniTn/DataMiningLab/exprMatTop{final_number_of_features}ByVar.tsv", sep='\t')
    print(matrix.shape)
    return matrix


if __name__ == '__main__':
    filter_samples = filter_generator(rf"C:/Users/aless/Documents/UniTn/DataMiningLab/objective.tsv")
    print(len(filter_samples))
    matrix = matrix_generator(barcodes_filter=filter_samples)
    print(matrix.shape)
    #matrix = matrix_trim(matrix, zero_ratio=0.9)
    #matrix = matrix_normalization(matrix, 'log')
    #print(matrix.iloc[1, :])
    #matrix = matrix.round(decimals=8)
    #print(matrix.iloc[1, :])
    #matrix.to_csv(r"C:/Users/aless/Desktop/expressionMatrix.csv")
    matrix = matrix_normalization(matrix, 'log')
    matrix = filter_by_variance(matrix)
    print(matrix.shape)

