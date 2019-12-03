import corex as ce
import pandas as pd
import matrix_generator as mg

N_FEATURES = mg.FINAL_NUMBER_OF_FEATURES
#expression_matrix_path = f'/home/alessio.albanese/exprMatTop{N_FEATURES}ByVar.tsv'
expression_matrix_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/exprMatTop{N_FEATURES}ByVar.tsv"
output_clusters_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/corex_clusters.tsv"
REPS = 3
N_CPU = 3
RAM = '2G'

expression_matrix = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, header=0)
expression_matrix = mg.filter_by_variance(expression_matrix, final_number_of_features=100)
print(expression_matrix.shape)
layer1 = ce.Corex(n_hidden=2, dim_hidden=1, marginal_description='gaussian', smooth_marginals=False,
                  n_repeat=REPS, n_cpu=N_CPU, ram=RAM).fit(expression_matrix)
clusters = layer1.clusters

