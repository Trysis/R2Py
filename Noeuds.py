# coding=utf-8
from ctypes import *
import math
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

# DataFrame 1 et 2 avec les sections nécessaires (ATOM/ HETATM)
ppdb_df1 = pd.concat([ppdb1.df[section] for section in records])
ppdb_df2 = pd.concat([ppdb2.df[section] for section in records])

# Dtf avec colonnes d'intérêts index, et columns_to_keep
ppdb_df1_to_keep = ppdb_df1[columns_to_keep].reset_index(level=0)[:500]
ppdb_df2_to_keep = ppdb_df2[columns_to_keep].reset_index(level=0)[:500]

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
NM_row = N * M

# Récupération dans un objet des fonctions dans C
f_node = CDLL('./edgex.so')

# Vertex
vertex = None
f_vertex = f_node.vertex
f_vertex.argtypes = [POINTER(POINTER(c_uint32)), c_char_p, c_char_p, c_size_t * 2, POINTER(c_size_t)]
f_vertex.restype = None

# Creation du vertex vide
vertex_xy = np.empty(shape=((NM_row), 2), dtype=np.uint32)

# Conversion en un type interpretable par ctypes
processV_xy = [POINTER(c_uint32)((c_uint32 * vertex_xy.shape[0])(*col_vertex)) for col_vertex in vertex_xy.T]
pointerV_xy = (POINTER(POINTER(c_uint32)))((POINTER(c_uint32) * vertex_xy.shape[1])(*processV_xy))

# Arguments à passer à notre fonction C
ppV_xy = pointerV_xy
X_elem_mem = c_char_p(str.encode("".join(X["element_symbol"])))
Y_elem_mem = c_char_p(str.encode("".join(Y["element_symbol"])))
XYsize_mem = (c_size_t * 2)(*[N, M])  # couple N, M spécifiant la taille de X et Y
vertex_size_mem = (c_size_t)(0)

# Lancement de la fonction sur C
f_vertex(ppV_xy, X_elem_mem, Y_elem_mem, XYsize_mem, byref(vertex_size_mem))

# On récupère le vertex
vertex_nrow = vertex_size_mem.value  # Nombre de ligne du vertex
vertex = np.array((ppV_xy[0][:vertex_nrow], ppV_xy[1][:vertex_nrow]), ndmin=2).T

# Edge
edge = None
f_edge = f_node.edge
f_edge.argtypes = [POINTER(POINTER(c_uint32)), POINTER(POINTER(c_uint32)), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)), c_size_t * 3, POINTER(c_size_t)]
f_edge.restype = None

# Nombre de lignes et colonnes du vertex
nrow_vertex = vertex.shape[0]
ncol_vertex = vertex.shape[1]

va, counts = np.unique(vertex.T[0], return_counts=True)

conca = np.concatenate([np.repeat(i, i) - np.arange(i) for i in counts])
NN_row_vertex = ((nrow_vertex - np.arange(nrow_vertex)).sum(dtype=np.uint64))  # taille de edge au départ
NN_row_vertex = int(NN_row_vertex)

if int(math.log10(NN_row_vertex) + 1) >= 7:
    NN_row_vertex = NN_row_vertex // 6

col_indices = X_dtf.columns.get_indexer(["x_coord", "y_coord", "z_coord"])

coordX = X_dtf.iloc[vertex.T[0], col_indices].T.values
coordY = Y_dtf.iloc[vertex.T[1], col_indices].T.values

# Creation de l'array edge vide de taille N*N (N = nrow(vertex))
# uint32 car max_atom_number = 99 999
edge_xy = np.empty(shape=((NN_row_vertex), 2), dtype=np.uint64)

# Conversion en un type interpretable par ctypes
processE_xy = [POINTER(c_uint32)((c_uint32 * edge_xy.shape[0])(*col_edge)) for col_edge in edge_xy.T]
pointerE_xy = POINTER(POINTER(c_uint32))((POINTER(c_uint32) * edge_xy.shape[1])(*processE_xy))

# Arguments à passer à notre fonction C
ppE_xy = pointerE_xy
ppV_xy

ncol_coord = 3

process_coordX = [cast((c_double * nrow_vertex)(*col_in_X), POINTER(c_double)) for col_in_X in coordX]
pointer_coordX = POINTER(POINTER(c_double))((POINTER(c_double) * ncol_coord)(*process_coordX))
pp_coordX = pointer_coordX #

process_coordY = [cast((c_double * nrow_vertex)(*col_in_Y), POINTER(c_double)) for col_in_Y in coordY]
pointer_coordY = POINTER(POINTER(c_double))((POINTER(c_double) * ncol_coord)(*process_coordY))
pp_coordY = pointer_coordY #

row_N_M_mem = (c_size_t * 3)(*[nrow_vertex, N, M])  # couple N, M spécifiant la taille de X et Y
edge_size_mem = (c_size_t)(NN_row_vertex)

f_edge(ppE_xy, ppV_xy, pp_coordX, pp_coordY, row_N_M_mem, byref(edge_size_mem))

edge_nrow = edge_size_mem.value  # Nombre de ligne du vertex
edge = np.array((ppE_xy[0][:edge_nrow], ppE_xy[1][:edge_nrow]), ndmin=2).T
print(edge)
print(f"{NN_row_vertex} {edge_nrow}")

exit()
