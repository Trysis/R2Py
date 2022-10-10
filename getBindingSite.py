# coding=utf-8
import sys
import os # gestion des fichiers
import argparse # gestion des arguments
from argparse import RawTextHelpFormatter

from Bio.PDB import PDBParser # Lecture fichier PDB
from Bio.PDB import NeighborSearch # Algorithme de NeighborSearch

# Gestion des arguments
# argparse permet de parser les options lors de l'exécution du script
# Nous exigereons d'exécuter notre fichier python en ajoutant en argument
#   obligatoire : le nom du fichier pdb à traiter
#   optionnelle : la localisation du répertoire en sortie (.../PDB_Atyping/)
# python pdb_file [-o] repository_out

args_parser = argparse.ArgumentParser(
    description="""
    Prend en argument un nom de fichier pdb, des chaînes cibles et le ligand
    Et renvoi les atomes (pour le moment les carbones alpha) qui sont à une
    distance inférieur à 5 Angström (paramètre modifiable par l'option -o)
    """,
    formatter_class=RawTextHelpFormatter
    )

args_parser.add_argument('pdb_file', type=str,help='Chemin vers le fichier pdb') # argument obligatoire

args_parser.add_argument('-o','--output', type=float, # argument optionnel
help="Chemin vers lequel nous enregistrons nos résultats")

args_parser.add_argument('-t','--tchain', type=str, # argument optionnel
help="Chaînes cibles pour l'évaluation des distances avec le ligand")

args_parser.add_argument('-l','--lchain', type=str, # argument optionnel
help="Chaîne du ligand dont on évaluera la distance avec nos chaînes cibles")

args_parser.add_argument('-c','--cutoff', type=float, # argument optionnel
help="Distance à laquelle on considère des atomes en contactes")

args = args_parser.parse_args() # Valeurs des arguments passés en argument par l'utilisateur
