# coding=utf-8

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

# Objets BioPandas
ppdb1 = PandasPdb()
ppdb2 = PandasPdb()

# Lecture des fichiers pdb
ppdb1.read_pdb('PDB/PDB_Binding/1AWR_P-A.pdb')
ppdb2.read_pdb('PDB/1BXB.pdb')

records = ['ATOM', 'HETATM'] # records à conserver dans le dataframe
columns_to_keep = ['element_symbol'] # Colonnes à conserver
columns_to_inner = ['element_symbol'] # Colonne à utiliser pour la jointure

# DataFrame 1 et 2 contenant les informations des sections section des
#   fichiers pdb
ppdb_df1 = pd.concat([ppdb1.df[section] for section in records])
ppdb_df2 = pd.concat([ppdb2.df[section] for section in records])

ppdb_df1_to_keep = ppdb_df1[columns_to_keep][:5] # Dtf avec colonnes d'intérêts
ppdb_df2_to_keep = ppdb_df2[columns_to_keep][:5] # Dtf avec colonnes d'intérêts

ppdb_df1_to_keep = ppdb_df1_to_keep.reset_index(level=0) # Ajout de la colonne indice
ppdb_df2_to_keep = ppdb_df2_to_keep.reset_index(level=0) # Ajout de la colonne indice

# Jointure entre les éléments du fichier PDB1 et PDB2
df_inner = ppdb_df1_to_keep.merge(ppdb_df2_to_keep, on=columns_to_inner, how='inner')

# Liste des indices présents dans la jointure entre PDB1 et PDB2
atm_index1 = df_inner["index_x"].unique() # Liste des indices chez PDB1
atm_index2 = df_inner["index_y"].unique() # Liste des indices chez PDB2

vertex = np.c_[df_inner["index_x"].T, df_inner["index_y"].T]
edge_list = []
edge = 0

def distance(atm_idx1, atm_idx2, ppdb_df1, ppdb_df2, coord_names = None):
    if coord_names == None:
        coord_names = ['x_coord', 'y_coord', 'z_coord']
    x = ppdb_df1[coord_names].iloc[[atm_idx1]].values
    y = ppdb_df2[coord_names].iloc[[atm_idx2]].values
    return(np.linalg.norm(x - y))

print(atm_index1)
print(atm_index2)
print(df_inner.head(20))