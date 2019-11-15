import numpy as np
import pandas as pd
from sklearn import manifold #''', decomposition'''
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

which_matrix = 'fpkms'
norm_method = 'log'
metadata_table = r"C:/Users/aless/Documents/DataMiningLab/clinical.tsv"


def log_norm(table):
    # log normalizes a table
    return np.log(table+0.00000001)


def fpkm_to_tpm_norm(table):
    # converts a (samples x genes) table of fpkms to tpms
    total_fpkms = list()
    for i, row in table.iterrows():
        tot_sample_fpkm = sum(row)
        total_fpkms.append(tot_sample_fpkm)
    table = table.div(total_fpkms, axis='rows')
    return table * (10 ** 6)


def deseq2_norm(table):
    # converts a row counts matrix to a DESeq2 normalized table
    mirror_table = np.log(np.float16(table.T))
    mirror_table = pd.DataFrame(data=mirror_table, index=table.T.index.values, columns=table.T.columns.values)
    any_zero_index = list()
    for i, row in mirror_table.iterrows():
        if row.isin([-np.inf]).any():
            any_zero_index.append(i)
        else:
            geo_mean = sum(row)/len(row)
            mirror_table.loc[i, :] = row-geo_mean
    mirror_table = mirror_table.drop(index=any_zero_index)
    scaling_factors = np.exp(mirror_table.median())
    return table.div(scaling_factors, axis='rows')


if which_matrix == 'raw_counts':
    matrix_path = r"C:/Users/aless/Documents/DataMiningLab/rawCountsExpressionMatrix.tsv"
elif which_matrix == 'fpkms':
    matrix_path = r"C:/Users/aless/Documents/DataMiningLab/expressionMatrix.tsv"
else:
    matrix_path = None
    exit()

samples = pd.read_csv(matrix_path, sep='\t', index_col=0, header=0)
samples = samples.loc[~(samples == 0).all(axis=1)]
samples = samples.T

# normalization of data:
if norm_method == 'deseq2':
    fig_name = "tsneDESeq2Norm.png"
    samples = deseq2_norm(samples)
elif norm_method == 'log':
    fig_name = "tsneLogNorm.png"
    samples = log_norm(samples)
elif norm_method == 'tpm':
    fig_name = "tsneTpmNorm.png"
    samples = fpkm_to_tpm_norm(samples)
elif norm_method == 'no_norm':
    fig_name = "tsneNoNorm.png"
else:
    print('Wrong Norm')
    exit()

samples = np.float32(samples)
(fig, subplots) = plt.subplots(1, 3, figsize=(18, 6))
perplexities = [40, 50, 60]
TSNEs = []
for i, perplexity in enumerate(perplexities):
    ax = subplots[i]
    tsne = manifold.TSNE(perplexity=perplexity, n_iter=10000, learning_rate=10.0)
    Y = tsne.fit(samples)
    print(f"Convergence with perplexity {perplexity} reached in {Y.n_iter_} iterations")
    Y = Y.embedding_
    TSNEs.append((Y,perplexity))

    ax.set_title("Perplexity=%d" % perplexity)
    ax.scatter(Y[:, 0], Y[:, 1])
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.axis('tight')
plt.savefig(fr"C:/Users/aless/Documents/DataMiningLab/{fig_name}")