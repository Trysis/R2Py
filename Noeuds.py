# coding=utf-8
from ctypes import *
import igraph as ig
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

# Objets BioPandas
ppdb1 = PandasPdb()
ppdb2 = PandasPdb()

# Lecture des fichiers pdb sous forme de DataFrame
ppdb1.read_pdb('PDB/2ptn.pdb')
ppdb2.read_pdb('PDB/1AWR.pdb')

records = ['ATOM', 'HETATM'] # records à conserver dans le dataframe
columns_to_keep = ["x_coord", "y_coord", "z_coord", 'element_symbol'] # Colonnes à conserver
columns_to_inner = ['element_symbol'] # Colonne à utiliser pour la jointure

# DataFrame 1 et 2 avec les sections nécessaires (ATOM/ HETATM)
ppdb_df1 = pd.concat([ppdb1.df[section] for section in records])
ppdb_df2 = pd.concat([ppdb2.df[section] for section in records])

# Dtf avec colonnes d'intérêts index, et columns_to_keep
ppdb_df1_to_keep = ppdb_df1[columns_to_keep].reset_index(level=0)
ppdb_df2_to_keep = ppdb_df2[columns_to_keep].reset_index(level=0)

# X sera le dataframe avec le plus grand nombre d'atomes
X, Y = None, None
X_dtf, Y_dtf = None, None

if ppdb_df1_to_keep.shape[0] >= ppdb_df2_to_keep.shape[0]:
    X = ppdb_df1_to_keep
    Y = ppdb_df2_to_keep
    X_dtf = ppdb_df1
    Y_dtf = ppdb_df2
else:
    X = ppdb_df2_to_keep
    Y = ppdb_df1_to_keep
    X_dtf = ppdb_df2
    Y_dtf = ppdb_df1

# Tailles N et M respectivement des dataframes X et Y
N = X.shape[0] # Nombre de lignes dans X
M = Y.shape[0] # Nombre de lignes dans Y

# Jointure entre les éléments du fichier PDB1 et PDB2
#   jointure sur les colonnes dans columns_to_inner
df_inner = X.merge(Y, on=columns_to_inner, how='inner')
I = df_inner.shape[0] # Nombre de lignes dans df_inner

X_idx = df_inner["index_x"]
Y_idx = df_inner["index_y"]

vertex = np.c_[X_idx, Y_idx] #?
distX = df_inner[["x_coord_x", "y_coord_x", "z_coord_x"]]
distY = df_inner[["x_coord_y", "y_coord_y", "z_coord_y"]]

func_edge = CDLL('D:/JM_Roude/Master_BioInformatique-Ingenieurie-de-Plateforme/UEs/Stage 2/R2Py/edgex.so')
func_edge.edge.restype = POINTER(POINTER(c_int))
# X et Y
n_col = (distX.shape[1])
n_row = (distX.shape[0])

pp_doubleX = POINTER(POINTER(c_double))
pp_doubleY = POINTER(POINTER(c_double))

process_distX = [cast((c_double * n_row)(*distX[col_in_X].values), POINTER(c_double)) for col_in_X in distX]
process_distY = [cast((c_double * n_row)(*distY[col_in_Y].values), POINTER(c_double)) for col_in_Y in distY]

X_pointer_val=(POINTER(c_double) * n_col)(*process_distX)
Y_pointer_val=(POINTER(c_double) * n_col)(*process_distY)

dist_pointerX = cast(X_pointer_val, pp_doubleX) #
dist_pointerY = cast(Y_pointer_val, pp_doubleY) #

X_idx_pointer = (c_int * n_row)(*X_idx.values)
Y_idx_pointer = (c_int * n_row)(*Y_idx.values)

nrow_edge = c_int()
edge_mat = POINTER(POINTER(c_int))
edge_mat.contents = func_edge.edge(byref(X_idx_pointer), byref(Y_idx_pointer),
                                    dist_pointerX, dist_pointerY,
                                    I, N, M, byref(nrow_edge))

edge = [(edge_mat.contents[0][i], edge_mat.contents[1][i]) for i in range(nrow_edge.value)]
print(edge)
exit()
