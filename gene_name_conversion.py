import pandas as pd
import mygene
mg = mygene.MyGeneInfo()

cluster_n = 0
clusters_file_path = rf"C:/Users/aless/Documents/UniTn/DataMiningLab/FinalReport/corex_clusters.tsv"
clusters = pd.read_csv(clusters_file_path, sep='\t', header=0, index_col=0)
cluster_one = clusters['clusters'] == cluster_n
print(clusters.loc[cluster_one, 'gene_ids'])
genes_list = list(map(lambda x: x.split('.')[0], clusters.loc[cluster_one, 'gene_ids']))
ens2sym = mg.querymany(genes_list, scopes='ensemblgene', fields='symbol', as_dataframe=True, returnall=False)
ens2sym = ens2sym.loc[ens2sym['notfound'] != True, 'symbol']
ens2sym.to_csv(rf"C:/Users/aless/Documents/UniTn/DataMiningLab/FinalReport/cluster{cluster_n+1}_sym.tsv", index=False)

