import corex as ce
import pandas as pd
import matrix_generator as mg

N_FEATURES = mg.FINAL_NUMBER_OF_FEATURES
expression_matrix_path = f'/home/alessio.albanese/exprMatTop{N_FEATURES}ByVar.tsv'
#expression_matrix_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/exprMatTop{N_FEATURES}ByVar.tsv"
output_clusters_path = f'/home/alessio.albanese/corex_clusters.tsv'
#output_clusters_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/corex_clusters.tsv"
output_custers_tcs_path = f'/home/alessio.albanese/corex_clusters_tcs.tsv'
#output_custers_tcs_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/corex_clusters_tcs.tsv"
REPS = 3
N_CPU = 100
RAM = 200

expression_matrix = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, header=0)
# expression_matrix = mg.filter_by_variance(expression_matrix, final_number_of_features=100)
layer1 = ce.Corex(n_hidden=150, dim_hidden=3, marginal_description='gaussian', smooth_marginals=True,
                  n_repeat=REPS, n_cpu=N_CPU, ram=RAM).fit(expression_matrix.T)
clusters = layer1.clusters
print(clusters)
print(len(clusters))
clusters_table = pd.DataFrame(list(zip(expression_matrix.index.values, clusters)), columns=['gene_ids', 'clusters'])
clusters_table.to_csv(output_clusters_path, sep='\t')
clusters_tcs = pd.DataFrame(list(zip(range(200), layer1.tcs)), columns=['clusters', 'tcs'])
clusters_tcs.to_csv(output_custers_tcs_path, sep='\t')
