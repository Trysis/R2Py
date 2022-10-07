# coding=utf-8

import sys, os, argparse
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB.SASA import ShrakeRupley

# Variables globales

ATOM_POS_RECORD = \
    {
    "ATOM":slice(0,6),"serial":slice(6,11),
    "atome_name":slice(12,16),"a_location":slice(16,17),
    "res_name":slice(17,20),"chain":slice(21,22),
    "res_nb":slice(22,26),"icode":slice(26,27),
    "x":slice(30,38),"y":slice(38,46),"z":slice(46,54),
    "occupancy":slice(54,60),"temperature":slice(60,66),
    "element":slice(76,78),"charge":slice(78,80)
    } # Position de chacune des colonnes dans un fichier pdb

# argparse module utile pour parser les options lors de l'exécution du script
# Nous exigereons d'exécuter notre fichier python en ajoutant en argument
#   le nom du fichier pdb à traiter

# python pdb_file [-o] repository_out
parser = argparse.ArgumentParser(description="")
parser.add_argument('pdb_file', type=str,help='Chemin vers le fichier pdb')
parser.add_argument('-o','--pdb_out', type=str,help="Emplacement où sera enregistré notre fichier de sortie")

args = parser.parse_args() # Valeurs des arguments passés en argument par l'utilisateur

# Fichier PDB
parser = PDBParser(QUIET=True)# Lecture Fichier PDB; QUIET=TRUE n'affiche pas msg d'erreur
io = PDBIO() # Sauvegarde fichier pdb seulements model/atomes
sr = ShrakeRupley() # Calcul de la surface accessible

# Vérification des chemins du fichier pdb, et répertoire

# PDB en entrée : Emplacement et nom de notre fichier pdb reçu en entrée
path_to_pdb = args.pdb_file.split("/")

pdb_directory = "./"
pdb_file = ""
pdb_file_name = ""
pdb_file_extension = ""

if len(path_to_pdb) > 1: # On assigne à pdb_directory le répertoire où se trouve notre fichier à changer
    pdb_directory = "/".join(path_to_pdb[:-1])+"/"

pdb_file_name, pdb_file_extension = os.path.splitext(path_to_pdb[-1]) #nom et extension du fichier pdb
pdb_file = path_to_pdb[-1] # Nom complet (avec extension) de notre fichier

# PDB en sortie : Emplacement et noms de notre répertoire et fichier pdb en sortie
pdb_out_directory = os.path.join(pdb_directory,"PDB_Atyping") # Chemin vers le répertoire à créer
pdb_out_file = pdb_file_name + "_Atyping" + pdb_file_extension # nom du fichier pdb à créer

# Chemin complets pour nos fichiers
path_file = args.pdb_file # chemin du fichier pdb en entrée
path_out_file = os.path.join(pdb_out_directory, pdb_out_file) # chemin du fichier pdb en sortie

# Création et vérification fichier/répertoire

#   Vérification de l'existence des fichier
path_exists = os.path.exists(path_file) # verifie si le fichier pdb existe
if not path_exists:
    error_message = \
    f"""
    Le fichier pdb dont vous avez spéficié le chemin n'éxiste pas
    path="{path_file}"\n
    """
    sys.exit(error_message)

#   Création de notre répertoire
if not os.path.isdir(pdb_out_directory): # verifie si le fichier pdb existe
    os.mkdir(pdb_out_directory)

# Biopython
# PDBParser
# Importation de notre fichier pdb dans un objet
structure = parser.get_structure(pdb_file_name, path_file) # PDBParser
io.set_structure(structure) # nous permettra de créer un fichier pdb

# Calcul de la surface accessible de chacun des atomes pour chaques chaînes de notre protéine
for model in structure:
    for chaine in model:
        sr.compute(chaine, level="A")

path_exists = os.path.exists(path_out_file) # verifie le chemin du fichier pdb modifie

exit()
# Fonctions

def is_Carbone_R(str_carbone): # Verifie si Carbone d'une chaine laterale
    str_carbone_c = str_carbone.strip()
    if len(str_carbone_c) <= 1: 
        #Verifie que nous ne sommes pas sur le carbone COO
        return False
    return "C" == str_carbone.strip()[0]

def is_Calpha(str_carbone):# verifie si Carbone alpha
    return "CA" == str_carbone.strip()

def is_Cbeta(str_carbone):# verifie si Carbone beta
    return "CB" == str_carbone.strip()

def is_aromatique(str_acide_amine): # Verifie si carbone aromatique
    aromatiques_valides = ["HIS","PHE","TYR","TRP"]
    return str_acide_amine.strip() in aromatiques_valides

#Ouverture du fichier pdb en lecture et ecriture pour ajouter les element ('a','b','A')
with open(path_out_file, "r+") as pdb_out:
    offset = 0 # position dans le cadre de lecture
    numero_ligne = 0 # numero de la ligne dans le fichier pdb

    # taille de la section element
    element_length = ATOM_POS_RECORD["element"].stop-ATOM_POS_RECORD["element"].start 
    
    for lines in pdb_out:
        if "ATOM" not in lines[ATOM_POS_RECORD["ATOM"]]: 
            pass

        atome_name = lines[ATOM_POS_RECORD["atome_name"]]# nom de l'atome a notre ligne
        res_name = lines[ATOM_POS_RECORD["res_name"]] # nom du residu a notre ligne

        if is_Calpha(atome_name):
            # On modifie le nom de l'élement à "a" si carbone alpha
            str_to_write = f"{'a':<{element_length}}"
            pdb_out.seek(offset+ATOM_POS_RECORD["element"].start+1)
            pdb_out.write(str_to_write)

        elif is_Cbeta(atome_name):
            # On modifie le nom de l'élement à "b" si carbone beta
            str_to_write = f"{'b':<{element_length}}"
            pdb_out.seek(offset+ATOM_POS_RECORD["element"].start+1)
            pdb_out.write(str_to_write)

        elif is_aromatique(res_name):
            if is_Carbone_R(atome_name):
                # Si residu aromatique, et carbone R, element devient "A"
                str_to_write = f"{'A':<{element_length}}"
                pdb_out.seek(offset+ATOM_POS_RECORD["element"].start+1)
                pdb_out.write(str_to_write)
        
        offset += len(lines) # Position de notre cadre de lecture
        pdb_out.seek(offset) # Necessaire si l'on veut maintenir la position de notre ligne

print("job done")