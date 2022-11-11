# coding=utf-8
from ctypes import *
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
columns_to_keep = ['element_symbol'] # Colonnes à conserver
columns_to_inner = ['element_symbol'] # Colonne à utiliser pour la jointure

# DataFrame 1 et 2 avec les sections nécessaires (ATOM/ HETATM)
ppdb_df1 = pd.concat([ppdb1.df[section] for section in records])[:6]
ppdb_df2 = pd.concat([ppdb2.df[section] for section in records])[:4]

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
distX = X_dtf.loc[X_idx, ["x_coord", "y_coord", "z_coord"]]
distY = Y_dtf.loc[Y_idx, ["x_coord", "y_coord", "z_coord"]]

func_edge = CDLL('D:/JM_Roude/Master_BioInformatique-Ingenieurie-de-Plateforme/UEs/Stage 2/R2Py/edgex.so')
func_edge.edge.restype = POINTER(POINTER(c_int))
# X et Y
n_col = (distX.shape[1])

pp_doubleX = POINTER(POINTER(c_double))
pp_doubleY = POINTER(POINTER(c_double))

process_distX = [cast((c_double * distX.shape[0])(*distX[col_in_X].values), POINTER(c_double)) for col_in_X in distX]
process_distY = [cast((c_double * distY.shape[0])(*distY[col_in_Y].values), POINTER(c_double)) for col_in_Y in distY]

X_pointer_val=(POINTER(c_double) * n_col)(*process_distX)
Y_pointer_val=(POINTER(c_double) * n_col)(*process_distY)

dist_pointerX = cast(X_pointer_val, pp_doubleX) #
dist_pointerY = cast(Y_pointer_val, pp_doubleY) #

X_idx_pointer = (c_int * distX.shape[0])(*X_idx.values)
Y_idx_pointer = (c_int * distY.shape[0])(*Y_idx.values)

print(distX)
print(distY)
nrow_edge = c_int()
edge_mat = POINTER(POINTER(c_int))
edge_mat.contents = func_edge.edge(byref(X_idx_pointer), byref(Y_idx_pointer),
                                    dist_pointerX, dist_pointerY,
                                    I, N, M, byref(nrow_edge))

print(nrow_edge.value)
print(f"{edge_mat.contents[0][30]} {edge_mat.contents[1][30]}")
edge = [(edge_mat.contents[0][i], edge_mat.contents[1][i]) for i in range(nrow_edge.value)]
exit()

# Edges
# N > M
# n inclus dans [0, N-1]; m inclus dans [0,M-1]
# i = 0
# for x_idx, y_idx in vertex:
#    dist_idx = distance_with_df(x_idx, y_idx, X_dtf, Y_dtf)
#    i += 1
#    for x_idx_prm, y_idx_prm in vertex[i:]:
#        dist_idx_prm = distance_with_df(x_idx_prm, y_idx_prm, X_dtf, Y_dtf)
#        md_idx = x_idx*N + y_idx
#        md_idx_prm = x_idx_prm*N + y_idx_prm
#
#        #print(f"(i={x_idx:.2f}, j={y_idx:.2f}) - (i'={x_idx_prm:.2f}, j'={y_idx_prm:.2f})'")
#        #print(f"md={md_idx:.2f} - md'={md_idx_prm:.2f}")
#
#        #print(f"i={md_idx//N} j={md_idx%N}")
#        #print(f"i'={md_idx_prm//N} j'={md_idx_prm%N}\n")
#        if(abs(dist_idx-dist_idx_prm) < seuil):
#            edge_list.append((md_idx, md_idx_prm))
#    
# print(edge_list)

# {Numéro indice: position dans le vertex}
atm_idx_pos =  {k: v[0] for k,v in X_idx.groupby(X_idx, sort=False).groups.items()}
atm_idx_keys = list(atm_idx_pos) #key

seuil = 0.1

edges = []

i = 0
j = 0
for x_idx, y_idx in vertex:
    key_idx = atm_idx_keys.index(x_idx) + 1
    if(key_idx >= len(atm_idx_keys)):
        break

    j = atm_idx_pos[atm_idx_keys[key_idx]]
    
    for x_idx_prm, y_idx_prm in vertex[j:]:
        dist_idx = distance_with_df(x_idx, x_idx_prm, X_dtf, X_dtf)
        dist_idx_prm = distance_with_df(y_idx, y_idx_prm, Y_dtf, Y_dtf)

        if(abs(dist_idx-dist_idx_prm) < seuil):
            edges.append((i, j))
        j += 1
    i += 1

