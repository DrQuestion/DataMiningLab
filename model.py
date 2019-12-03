import pandas
import numpy as np

"""
    Select a subset of table with certain genes
    parameters:
        table: pandas dataframe with rows are samples, columns are genes
        sets: list of sets in which each set contains interested genes
"""
def select_table(table, sets):
    union_set = {}
    for s in sets:
        union_set = union_set.union(s)
    return table[list(union_set)]
