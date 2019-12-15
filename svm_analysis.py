import numpy as np
import pandas as pd
import model
import nessra_genes_sets_generator as ng
import mygene 


#------------- MATRIX and CLASS VECTOR
matrix = pd.read_csv("~/Desktop/expressionMatrix.tsv", sep = "\t", header=0)
class_vector = pd.read_csv("~/GoogleDrive/DataMiningLab/Model/objective.tsv", sep = "\t", header=0)
print("Matrix shape: {}".format(matrix.shape))
print(matrix.head(2))
print("")

# transposing matrix
matrix_t = matrix.T
matrix_t.columns = matrix_t.iloc[0]
matrix_t = matrix_t.drop("Gene_ID", axis=0)

# Resetting column names without the version of ENS
col_names = list(matrix_t.columns.values)
new_col_names = []
for name in col_names:
    temp_split = name.split(".")
    new_col_names.append(temp_split[0])
matrix_t.columns = new_col_names
print("Transpose matrix shape: {}".format(matrix_t.shape))
print("")

# class vector 
print("Objective shape: {}".format(class_vector.shape))    

# generating barcode of samples
samples_bar = []
for name in list(matrix_t.index.values):
    samples_bar.append(name[:-4])
fil_for_class = class_vector.barcode.isin(samples_bar)
class_vector = class_vector[fil_for_class]
print("Objective shape after filter for matrix samples: {}".format(class_vector.shape))

# Final vector for classifying
train = class_vector["over2years"]
print("Final vector for classiftibg is a {}, first elements: {}".format(type(train), train[:5]))
print("")


#----------------- NESSRA -----------------
print("----SVM on NESSRA----")
print("")
# importing files output
list_apc = [pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160679_Hs.expansion', sep = ",", header=1),
            pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160679_Hs.interactions', sep = ",", header=1)]
list_smad4 = [pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160680_Hs.expansion', sep = ",", header=1),
            pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160680_Hs.interactions', sep = ",", header=1)]
list_kras = [pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160812_Hs.expansion', sep = ",", header=1),
            [None]]
list_braf = [pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/147717_Hs.expansion', sep = ",", header=1),
            [None]]

# importing tcode conversion matrix
tcode_mat = pd.read_csv('~/GoogleDrive/DataMiningLab/hgnc_filtered_anno.csv', sep = ",", header=0)

# Set filter tresholds and max genes per set here!
freq_treshold = 0.9
n_genes_per_file = 450
tcode_symbol_dict = ng.tcodeSymbolDictGen(tcode_mat)
set_apc = ng.setGen(list_apc[0], list_apc[1], freq_treshold, n_genes_per_file, tcode_symbol_dict, useIntFile=True)
set_smad4 = ng.setGen(list_smad4[0], list_smad4[1], freq_treshold, n_genes_per_file, tcode_symbol_dict, useIntFile=False)
set_kras = ng.setGen(list_kras[0], list_kras[1], freq_treshold, n_genes_per_file, tcode_symbol_dict, useIntFile=False)
set_braf = ng.setGen(list_braf[0], list_braf[1], freq_treshold, n_genes_per_file, tcode_symbol_dict, useIntFile=False)

list_of_sets = [set_apc, set_smad4, set_braf, set_braf ]

print("Len of set apc: {}".format(len(set_apc)))
print("Len of set smad4: {}".format(len(set_smad4)))
print("Len of set kras: {}".format(len(set_kras)))
print("Len of set braf: {}".format(len(set_braf)))
print("Subsetting matrix with NESSRA output genes")
ns_mat = model.select_table(matrix_t, list_of_sets)
print("Transpose matrix filtered shape: {}".format(ns_mat.shape))
print(ns_mat.head(2))
print("")

#------------ SVM ANALYSIS ------------
# Import train_test_split function
from sklearn.model_selection import train_test_split

# Split dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(ns_mat, train, test_size=0.2,random_state=300) # 80% training and 20% test

# train the model
clf_ns = model.training(X_train, y_train)

y_pred = clf_ns.predict(X_test)
#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics

print("")
print("--- Evaluating test set ---")
print("Test matrix shape: {} \n{} \netc. \nTest class vector: {} \n{}, etc.\n".format(X_test.shape, X_test.head(2), len(y_test), y_test[:20]))
print("Frequency filter used in NESSRA output: {}, final number of genes: {}". format(freq_treshold, X_test.shape[1]))
print("")
# Model Accuracy: how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

# Model Precision: what percentage of positive tuples are labeled as such?
print("Precision:",metrics.precision_score(y_test, y_pred))

# Model Recall: what percentage of positive tuples are labelled as such?
print("Recall:",metrics.recall_score(y_test, y_pred))


# Save genes as SYMBOL

ens_gene_list = list(ns_mat.columns.values)
mg = mygene.MyGeneInfo()
symbol_genes = mg.querymany(ens_gene_list, scopes='ensembl.gene', fields='symbol', species='human')
gene_list = []
for el in symbol_genes:
    if "symbol" in el:
        gene_list.append(el["symbol"])
print("len of ens genes: {}, len of symbol genes: {}".format(len(ens_gene_list), len(gene_list)))
print(', '.join(gene_list))

