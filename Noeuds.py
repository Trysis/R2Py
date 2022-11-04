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
M = df_inner.shape[0] # Nombre de lignes dans df_inner
X_idx = df_inner["index_x"]
Y_idx = df_inner["index_y"]

vertex = np.c_[X_idx, Y_idx] #?

func_edge = CDLL('edgex.so')
# func_edge.edge.argtypes = POINTER(POINTER(c_int)), POINTER(c_int)
# func_edge.edge.restype = None

val = [(c_int * M)(*X_idx.values), (c_int * M)(*Y_idx.values)]
a = POINTER(c_int)
a.contents = (c_int * M)(*X_idx.values)

b = POINTER(c_int)
b.contents = (c_int * M)(*Y_idx.values)

mem = POINTER(POINTER(c_int))()
mem.contents = c_int * 2
mem[0].contents = a
mem[1].contents = b

func_edge.edge(mem, M)

exit()

# Calcul la distance entre 2 atomes à partir des
def distance_with_df(atm_idx1, atm_idx2, ppdb_df1, ppdb_df2):
    """
    Calcul la distance entre 2 atomes à partir
      des coordonnées x, y, z issus des dataframes

    Arguments :
        atm_idx1 : indice de l'atome 1
        atm_idx2 : indice de l'atome2
        ppdb_df1 : dataframe 1
        ppdb_df2 : dataframe 2
    """
    
    coord_names = ['x_coord', 'y_coord', 'z_coord']
    x = ppdb_df1[coord_names].iloc[[atm_idx1]].values
    y = ppdb_df2[coord_names].iloc[[atm_idx2]].values
    return(np.linalg.norm(x - y))


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

