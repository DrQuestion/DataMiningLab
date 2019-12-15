import textwrap
import corex as ce
import vis_corex as vc
import pandas as pd
import matrix_generator as mg

N_FEATURES = mg.FINAL_NUMBER_OF_FEATURES

N_HIDDEN = 200
REPS = 3
N_CPU = 20
RAM = 30

max_edges = 100

expression_matrix_path = f'/home/alessio.albanese/exprMatTop{N_FEATURES}ByVar.tsv'
#expression_matrix_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/exprMatTop{N_FEATURES}ByVar.tsv"
output_clusters_path = f'/home/alessio.albanese/corex_clusters.tsv'
#output_clusters_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/corex_clusters.tsv"
output_custers_tcs_path = f'/home/alessio.albanese/corex_clusters_tcs.tsv'
#output_custers_tcs_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/corex_clusters_tcs.tsv"
output_graph_path = '/home/alessio.albanese/graph_prune_' + str(max_edges)
#output_graph_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/corex_graph" + str(max_edges)
output_mis_path = rf"/home/alessio.albanese/corex_mis.tsv"
#output_mis_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/corex_mis.tsv"

expression_matrix = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, header=0)
# expression_matrix = mg.filter_by_variance(expression_matrix, final_number_of_features=100)
layer1 = ce.Corex(n_hidden=N_HIDDEN, dim_hidden=3, marginal_description='gaussian', smooth_marginals=True,
                  n_repeat=REPS, n_cpu=N_CPU, ram=RAM).fit(expression_matrix.T)
layer2 = ce.Corex(n_hidden=1, dim_hidden=3, marginal_description='discrete', smooth_marginals=True,
                  n_repeat=REPS, n_cpu=N_CPU, ram=RAM).fit(layer1.labels)
corexes = [layer1, layer2]

clusters = layer1.clusters
clusters_table = pd.DataFrame(list(zip(expression_matrix.index.values, clusters)), columns=['gene_ids', 'clusters'])
clusters_table.to_csv(output_clusters_path, sep='\t')
clusters_tcs = pd.DataFrame(list(zip(range(N_HIDDEN), layer1.tcs)), columns=['clusters', 'tcs'])
clusters_tcs.to_csv(output_custers_tcs_path, sep='\t')
mis = layer1.mis
mis_table = pd.DataFrame(mis, columns=expression_matrix.index.values)
mis_table.to_csv(output_mis_path, sep='\t')


# Build the graph representation file .dot
column_label = expression_matrix.index.values
column_label = list(map(lambda q: '\n'.join(textwrap.wrap(q, width=17, break_long_words=False)), column_label))
weights = [corex.alpha[:, :, 0].clip(0, 1) * corex.mis for corex in corexes]
node_weights = [corex.tcs for corex in corexes]
g = vc.make_graph(weights, node_weights, column_label, max_edges=max_edges)
h = g.copy()
vc.edge2pdf(h, output_graph_path, labels='label', directed=True, makepdf=False)



