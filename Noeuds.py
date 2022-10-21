import numpy as np
import itertools
from itertools import accumulate, repeat, chain
from itertools import product
from Bio.PDB import PDBParser # Lecture fichier PDB



# From https://stackoverflow.com/questions/11144513/
# Produit cartésien à N dimensions

def cartesian_product_pp(arrays):
    la = len(arrays)
    L = *map(len, arrays), la # Tailles de nos arrays et nb arrays
    dtype = np.result_type(*arrays)
    arr = np.empty(L, dtype=dtype) # Creation de notre array
    arrs = *accumulate(chain((arr,), repeat(0, la-1)), np.ndarray.__getitem__),
    idx = slice(None), *itertools.repeat(None, la-1)
    for i in range(la-1, 0, -1):
        arrs[i][..., i] = arrays[i][idx[:la-i]]
        arrs[i-1][1:] = arrs[i]
    arr[..., 0] = arrays[0][idx]
    return arr.reshape(-1, la)

def cartesian_product_same_elm(atm_iter1, atm_iter2):
    return [(x,y) for x, y in product(atm_iter1, atm_iter2) if x.element == y.element]

def cartesian_product_same_elm4(atm_iter1, atm_iter2):
    N, M, e = len(atm_iter1), len(atm_iter2), 0
    V = np.empty((N*M,2), dtype=object)
    for i in atm_iter2:
        for j in atm_iter1:
            if(i.element == j.element):
                V[e,] = i,j
                e += 1
    return V[:e,]

def cartesian_product_same_elm5(atm_iter1, atm_iter2):
    N, M, e = len(atm_iter1), len(atm_iter2), 0
    V = np.empty((N*M,2))
    for i in range(M):
        for j in range(N):
            if(atm_iter2[i].element == atm_iter1[j].element):
                V[e,] = i+1,j+1
                e += 1
    return V[:e,]

def cartesian_product_same_elm_withSorted(atm_iter1, atm_iter2):
    N, M, e = len(atm_iter1), len(atm_iter2), 0
    V = np.empty((N*M,2))
    j = 0

    for i in range(M):
        while(j < N):
            if(atm_iter2[i].element != atm_iter1[j].element):
                
                break
            V[e,] = i+1,j+1
            e += 1
            j += 1

    return V[:e,]


parser = PDBParser(QUIET=True) # Lecture Fichier PDB; QUIET=TRUE n'affiche pas de msg d'erreur
structure = parser.get_structure("nimp", "PDB/1BXB.pdb") # PDBParser
model0 = structure[0]
atoms = sorted([i for i in model0.get_atoms()][:10000], key=lambda x:x.element)
atoms_2 = atoms.copy()

prod = cartesian_product_pp((atoms,atoms_2))
