# coding=utf-8

import sys
import os # gestion des fichiers
import argparse # gestion des arguments
from argparse import RawTextHelpFormatter

from Bio.PDB import Selection
from Bio.PDB import PDBParser # Lecture fichier PDB
from Bio.PDB import NeighborSearch # Algorithme de NeighborSearch

# Gestion des arguments
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

args_parser.add_argument('-c','--cutoff', type=float, # argument optionnel,
default=5, help="Distance à laquelle on considère des atomes en contactes")

# Récupération des arguments
args = args_parser.parse_args() # Valeurs des arguments passés en argument par l'utilisateur

path_file = args.pdb_file # Nom du chemin vers le fichier pdb
out_args = args.output # Chemin de destination du fichier en sortie

# Vérification de l'existence des fichier
path_exists = os.path.exists(path_file) # verifie si le fichier pdb existe
if not path_exists:
    error_message = \
    f"""
    Le fichier pdb dont vous avez spéficié le chemin n'éxiste pas
    path="{path_file}"\n
    """
    sys.exit(error_message)


# Vérification des chemins du fichier pdb, et répertoire
# PDB en entrée : Emplacement et nom de notre fichier pdb reçu en entrée
path_to_pdb = path_file.split("/") # Liste des noms séparés par le slash

pdb_directory = "./" # Correspondra à l'emplacement du répertoire du fichier pdb
pdb_file = "" # nom du fichier pdb complet
pdb_file_name = "" # nom du fichier (sans extension)
pdb_file_extension = "" # extension du fichier

if len(path_to_pdb) > 1:
    # On assigne à pdb_directory le répertoire où se trouve notre fichier
    pdb_directory = "/".join(path_to_pdb[:-1])+"/"

pdb_file = path_to_pdb[-1] # Nom complet (avec extension) de notre fichier
pdb_file_name, pdb_file_extension = os.path.splitext(path_to_pdb[-1]) # nom et extension du fichier pdb

# PDB en sortie : Emplacement et noms de notre répertoire et fichier pdb en sortie
pdb_out_directory_name = "PDB_Binding" # Nom du répertoire que l'on va créer
pdb_out_directory = os.path.join(pdb_directory,pdb_out_directory_name)

if out_args != None: # | Chemin vers le répertoire à créer
    pdb_out_directory = os.path.join(out_args,pdb_out_directory_name)

# pdb_out_file, nom du fichier à voir
pdb_out_file = pdb_file_name + "_Binding" + pdb_file_extension # nom du fichier pdb à créer

# Chemin complets pour nos fichiers
path_out_file = os.path.join(pdb_out_directory, pdb_out_file) # chemin du fichier pdb en sortie

# Création de notre répertoire pour la conservation du fichier pdb sortant
#if not os.path.isdir(pdb_out_directory): # verifie si le répertoire n'existe pas déjà
    #os.mkdir(pdb_out_directory)

tchains_args = args.tchain # Chaînes cible 
ligands_args = args.lchain # Ligand cible
cutoff_args = args.cutoff # Distance à laquelle on considère un contact
# seuil pour cutoff ?

tchains_liste = tchains_args.replace(","," ").split()
ligands_liste = ligands_args.replace(","," ").split()

# Biopython
# PDBParser
# Lecture des données PDB à l'aide de PDBParser

parser = PDBParser(QUIET=True) # Lecture Fichier PDB; QUIET=TRUE n'affiche pas de msg d'erreur
structure = parser.get_structure(pdb_file_name, path_file) # PDBParser

# On vérifie dans un premier temps que les chaînes passés en arguments
# sont des chaînes qui existent dans notre fichier pdb
chaines_id_in = set(tchains_liste + ligands_liste) # chaines passés en arguments
chaines_id_pdb = {chaine_id.id for chaine_id in structure[0].get_chains()} # chaines du fichier pdb
chaines_diff = chaines_id_in.difference(chaines_id_pdb) # Chaines non existantes dans notre fichier pdb

# Verification chaînes identiques ?
# Si l'utilisateur a spécifié des chaînes non existantes
if len(chaines_diff)>0:
    error_message = \
        f"""
        Le(s) chaînes {chaines_diff} n'existe(nt) pas\n
        """
    sys.exit(error_message)

dictionnaire = {}
atomes_du_ligand = [x for x in structure[0]["A"].get_atoms()]
atomes_cibles = [x for x in structure[0]["A"].get_atoms()]

neighbor = NeighborSearch(atomes_du_ligand)
print(neighbor.search(atomes_cibles[0].coord,cutoff_args,'A'))

print(f"chaines = {tchains_liste}")
print(f"ligands = {ligands_liste}")
print(f"cutoff = {cutoff_args}")