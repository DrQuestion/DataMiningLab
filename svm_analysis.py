import numpy as np
import pandas as pd
import model
import nessra_genes_sets_generator as ng


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
print(matrix_t.head(2))
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
train = class_vector.iloc[:, -1].tolist()
print("Final vector for classiftibg is a {}, first elements: {}".format(type(train), train[:5]))
print("")


#----------------- NESSRA -----------------
print("----SVM on NESSRA----")
print("")
# importing files output
list_79 = [pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160679_Hs.expansion', sep = ",", header=1),
            pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160679_Hs.interactions', sep = ",", header=1)]
list_80 = [pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160680_Hs.expansion', sep = ",", header=1),
            pd.read_csv('~/GoogleDrive/DataMiningLab/Nessra/160680_Hs.interactions', sep = ",", header=1)]
# importing tcode conversion matrix
tcode_mat = pd.read_csv('~/GoogleDrive/DataMiningLab/hgnc_filtered_anno.csv', sep = ",", header=0)

# Set filter tresholds and max genes per set here!
freq_treshold = 0.3
n_genes_per_file = 50
tcode_symbol_dict = ng.tcodeSymbolDictGen(tcode_mat)
set79 = ng.setGen(list_79[0], list_79[1], freq_treshold, n_genes_per_file, tcode_symbol_dict, useIntFile=False)
set80 = ng.setGen(list_80[0], list_80[1], freq_treshold, n_genes_per_file, tcode_symbol_dict, useIntFile=False)
list_of_sets = [set79, set80]

print("Subsetting matrix with NESSRA output genes")
ns_mat = model.select_table(matrix_t, list_of_sets)
print("Transpose matrix filtered shape: {}".format(ns_mat.shape))
print(ns_mat.head(2))
print("")

#------------ SVM ANALYSIS ------------
# Import train_test_split function
from sklearn.model_selection import train_test_split

# Split dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(ns_mat, train, test_size=0.2,random_state=109) # 80% training and 20% test

# train the model

clf_ns = model.training(X_train, y_train)

y_pred = clf_ns.predict(X_test)
#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics

print("")
print("--- Evaluating test set ---")
print("Test matrix shape: {} \n{} \netc. \nTest class vector: {} \n{}, etc.\n".format(X_test.shape, X_test.head(2), len(y_test), y_test[:20]))
# Model Accuracy: how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

# Model Precision: what percentage of positive tuples are labeled as such?
print("Precision:",metrics.precision_score(y_test, y_pred))

# Model Recall: what percentage of positive tuples are labelled as such?
print("Recall:",metrics.recall_score(y_test, y_pred))
