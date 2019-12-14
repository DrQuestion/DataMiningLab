import numpy as np
import pandas as pd
import mygene



def tcodeSymbolDictGen(tcode_mat):
    '''
    Takes tcode matrix and creates a dictionary of tcode to ENS
    '''
    tcode_symbol_dict = {}
    print("I'm importing tcode-to-Symbol dictionary, this will take a minute")
    for i, row in tcode_mat.iterrows():
        tcode_symbol_dict[row[0].lower()] = row[10]
    print("Done, dictionary ready")
    
    return tcode_symbol_dict

def tcodeToEns(t_gene, tcodeSymbol_dict):
    '''
    Takes tcode_gene and dictionary tcode to symbol, returns ens using package mygene!
    '''
    symbol_gene = tcodeSymbol_dict[t_gene]

    return symbol_gene


def setGen(exp_file, int_file, freq_tres, n_genes, tcode_symbol_dict):
    '''
    Takes as input expansion ordered list of genes, related interaction file, frequency treshold for filter, max number of genes per file
    Returns a set of genes
    '''
  
    symbol_genes = []
    # Take genes from expansion file with filter treshold for frequency and converts them into Symbol
    for i, gene in exp_file.iterrows():
        if len(symbol_genes) < n_genes:
            symbol_gene = tcodeToEns(gene[1],  tcode_symbol_dict)
            if symbol_gene not in symbol_genes and gene[3] > freq_tres:
                symbol_genes.append(symbol_gene)
        else:
            break    
    
    # Check for final length, if not n_genes look into interactions
    if len(symbol_genes) < n_genes:
        for i, ints in int_file.iterrows():
            if len(symbol_genes) < n_genes and ints[4] > freq_tres:
                # check first x position
                gen_posX = tcodeToEns(ints[1],  tcode_symbol_dict)
                gen_posY = tcodeToEns(ints[2],  tcode_symbol_dict)
                if gen_posX in symbol_genes:
                    if gen_posY not in symbol_genes:
                        symbol_genes.append(gen_posY)
                # eventually check y position
                elif gen_posY in symbol_genes:
                    if gen_posX not in symbol_genes:
                        symbol_genes.append(gen_posX)
            else:
                break

    # Convert symbol genes into ens and add them into set, return set
    gene_set = set()
    mg = mygene.MyGeneInfo()
    ens_genes = mg.querymany(symbol_genes, scopes='symbol', fields='ensembl.gene', species='human')
    for el in ens_genes:
        if "ensembl" in el:
            if len(el["ensembl"]) > 1:
                for gene in el["ensembl"]:
                    gene_set.add(gene["gene"])
            else:
                gene_set.add(el["ensembl"]["gene"])
        
    return gene_set

##### Running functions #####


if __name__ == '__main__':

    list_79 = [pd.read_csv('~/GoogleDrive/DataMining/Nessra/160679_Hs.expansion', sep = ",", header=1),
            pd.read_csv('~/GoogleDrive/DataMining/Nessra/160679_Hs.interactions', sep = ",", header=1)]
    list_80 = [pd.read_csv('~/GoogleDrive/DataMining/Nessra/160680_Hs.expansion', sep = ",", header=1),
            pd.read_csv('~/GoogleDrive/DataMining/Nessra/160680_Hs.interactions', sep = ",", header=1)]
    tcode_mat = pd.read_csv('~/GoogleDrive/DataMining/hgnc_filtered_anno.csv', sep = ",", header=1)

    # Set filter tresholds and max genes per set here!
    freq_treshold = 0.3
    n_genes_per_file = 40
    tcode_symbol_dict = tcodeSymbolDictGen(tcode_mat)

    set79 = setGen(list_79[0], list_79[1], freq_treshold, n_genes_per_file, tcode_symbol_dict)
    set80 = setGen(list_80[0], list_80[1], freq_treshold, n_genes_per_file, tcode_symbol_dict)

    list_of_sets = [set79, set80]






