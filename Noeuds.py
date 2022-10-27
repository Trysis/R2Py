# coding=utf-8

import numpy as np
from biopandas.pdb import PandasPdb

ppdb = PandasPdb()
ppdb.read_pdb('PDB/1BXB.pdb')

ppdb2 = PandasPdb()
ppdb2.read_pdb('PDB/1BXB.pdb')

columns_k = ['atom_number','element_symbol']
columns_inner = ['element_symbol']

ppdb_df1 = ppdb.df['ATOM'][columns_k][:5]
ppdb_df2 = ppdb2.df['ATOM'][columns_k][:10]

mergedStuff = ppdb_df1.merge(ppdb_df2, on=columns_inner, how='inner')
a = mergedStuff["atom_number_x"].unique()