# coding=utf-8

import sys
import os # gestion des fichiers
import argparse # gestion des arguments
from argparse import RawTextHelpFormatter
from Bio.PDB import PDBParser # Lecture fichier PDB
from Bio.PDB import NeighborSearch # Algorithme de NeighborSearch

# Gestion des arguments passés par l'utilisateur
args_parser = argparse.ArgumentParser(
    description="""
    Prend en argument un nom de fichier pdb, des chaînes cibles et le ligand
    Et renvoi les atomes (pour le moment les carbones alpha) qui sont à une
    distance inférieur à 5 Angström de notre ligand
    """,
    formatter_class=RawTextHelpFormatter
    )

args_parser.add_argument('pdb_file', type=str,help='Chemin vers le fichier pdb') # argument obligatoire

args_parser.add_argument('-o','--output', type=float, # argument optionnel
help="Chemin vers lequel nous enregistrons nos résultats")

args_parser.add_argument('-t','--tchain', type=str, # argument optionnel
help="Chaînes cibles pour l'évaluation des distances avec le ligand")

args_parser.add_argument('-l','--lchain', type=str, # argument "optionnel"
help="Chaîne du ligand dont on évaluera la distance avec nos chaînes cibles")

args_parser.add_argument('-c','--cutoff', type=float, # argument optionnel,
default=5, help="Distance en Angström à laquelle on considère"
            +" des atomes en contactes, vaut 5 par défaut")

# Valeurs des arguments passés par l'utilisateur
args = args_parser.parse_args() # Récupération des arguments

path_in = args.pdb_file # Nom du chemin vers le fichier pdb
path_out = args.output # Chemin de destination du fichier en sortie
tchains = args.tchain # Chaînes cible 
lchains = args.lchain # Ligand cible
cutoff = args.cutoff # Distance à laquelle on considère un contact

# Vérification de l'existence du fichier pdb
path_exists = os.path.exists(path_in) # renvoi True si le fichier existe
if not path_exists:
    error_message = \
    f"""
    Le fichier pdb dont vous avez spéficié le chemin n'éxiste pas
    path="{path_in}"\n
    """
    sys.exit(error_message)

# Vérifie si l'utilisateur a passé en argument une valeur
# pour la chaîne du ligand
if lchains == None: 
    error_message = \
    f"""
    Veuillez spécifier la chaîne du ligand option [-l]
    """
    sys.exit(error_message)

# PDB en entrée : Emplacement et nom de notre fichier pdb reçu en entrée
path_to_pdb = path_in.split("/") # Liste des noms séparés par le slash

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

if path_out != None: # | Chemin vers le répertoire à créer
    pdb_out_directory = os.path.join(path_out,pdb_out_directory_name)

# seuil pour cutoff ?


# Biopython
# PDBParser
# Lecture des données PDB à l'aide de PDBParser

parser = PDBParser(QUIET=True) # Lecture Fichier PDB; QUIET=TRUE n'affiche pas de msg d'erreur
structure = parser.get_structure(pdb_file_name, path_in) # PDBParser

# Récupération des chaînes passées en argument
tchains_set = set(tchains.replace(","," ").split()) # Target chains
lchains_set = set(lchains.replace(","," ").split()) # Ligand chains

# On vérifie dans un premier temps que les chaînes passées en arguments
# sont des chaînes qui existent dans notre fichier pdb
chaines_id_in = tchains_set | ligands_liste # Union des chaines chez tchains et lchains
chaines_id_pdb = {chaine_id.id for chaine_id in model for model in structure} # chaines du fichier pdb
chaines_diff = chaines_id_in.difference(chaines_id_pdb) # Chaines non existantes dans notre fichier pdb

# Verification chaînes identiques ?
# Si l'utilisateur a spécifié des chaînes non existantes
if len(chaines_diff)>0:
    error_message = \
        f"""
        Le(s) chaînes {chaines_diff} n'existe(nt) pas\n
        """
    sys.exit(error_message)

# Création de notre répertoire pour la conservation du fichier pdb sortant
#if not os.path.isdir(pdb_out_directory): # verifie si le répertoire n'existe pas déjà
    #os.mkdir(pdb_out_directory)
# pdb_out_file, (nom du fichier à voir)
pdb_out_file = pdb_file_name + "_Binding" + pdb_file_extension # nom du fichier pdb à créer

# Chemin complets pour nos fichiers
path_out_file = os.path.join(pdb_out_directory, pdb_out_file) # chemin du fichier pdb en sortie

dictionnaire = {}
atomes_du_ligand = [x for x in structure[0]["A"].get_atoms()]
atomes_cibles = [x for x in structure[0]["A"].get_atoms()]

neighbor = NeighborSearch(atomes_du_ligand)
print(neighbor.search(atomes_cibles[0].coord,cutoff,'A'))

print(f"chaines = {tchains_liste}")
print(f"ligands = {ligands_liste}")
print(f"cutoff = {cutoff}")