import os
import re
import numpy as np
import pandas as pd

def setGen(exp_file, int_file, freq_tres, n_genes):
    '''
    Takes as input expansion ordered list of genes, related interaction file, frequency treshold for filter, max number of genes per file
    Returns a set of genes
    '''
    gene_set = set()

    # Take genes from expansion file with filter treshold for frequency
    for i, gene in exp_file.iterrows():
        if len(gene_set) < n_genes:
            if gene[1] not in gene_set and gene[3] > freq_tres:
                gene_set.add(gene[1])
        else:
            break    

    # Check for final length, if not n_genes look into interactions
    if len(gene_set) < n_genes:
        for i, ints in int_file.iterrows():
            if len(gene_set) < n_genes and ints[4] > freq_tres:
                # check first x position
                if ints[1] in gene_set:
                    if ints[2] not in gene_set:
                        gene_set.add(ints[2])
                # eventually check y position
                elif ints[2] in gene_set:
                    if ints[1] not in gene_set:
                        gene_set.add(ints[1])
            else:
                break
    
    return gene_set


# Imports all file from working directory and divides them into expansion and interaction files
folder_files = os.listdir('./')
exp_files = [name for name in folder_files if re.search('expansion', name)]
int_files = [name for name in folder_files if re.search('interaction', name)]


# Set filter tresholds and max genes per set here!
freq_treshold = 0.3
n_genes_per_file = 40

list_of_sets = []
for i in range(len(exp_files)):
    exp_genes = pd.read_csv(exp_files[i], sep = ",", header=1)
    int_genes = pd.read_csv(int_files[i], sep = ",", header=1)
    # apply function
    my_set = setGen(exp_genes, int_genes, freq_treshold, n_genes_per_file)
    print("looking at files {} and {}".format(exp_files[i], int_files[i]))
    print("len of gene set: {}".format( len(my_set)))
    print("")
    list_of_sets.append(my_set)






