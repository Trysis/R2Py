# coding=utf-8

import sys
import os # gestion des fichiers
import argparse # gestion des arguments
from argparse import RawTextHelpFormatter
from Bio.PDB import PDBParser # Lecture fichier PDB
from Bio.PDB.SASA import ShrakeRupley # Algorithme de calcul de la surface accessible

# Import nos fichiers
import auxiliaires as aux

### Fonctions
# Vérification de l'atome
def is_Carbone(atome_str):
    """Renvoi True si atome_str est un carbone"""
    return "C" == atome_str[0]


def is_Azote(atome_str):
    """Renvoi True si atome_str est un azote"""
    return "N" == atome_str


def is_Oxygen(atome_str):
    """Renvoi True si atome_str est un oxygène"""
    return "O" == atome_str


def is_Calpha(atome_str):
    """Renvoi True si atome_str est un carbone alpha"""
    return "CA" == atome_str


def is_Cbeta(atome_str):
    """Renvoi True si atome_str est un carbone bétâ"""
    return "CB" == atome_str


def is_CAromatique(str_acide_amine, atome_str):
    """Renvoi True si atome_str est un carbone est aromatique"""
    aromatiques = {
        "HIS":["CG","CD2","CE1"],
        "PHE":["CG","CD1","CD2","CE1","CE2","CZ"],
        "TYR":["CG","CD1","CD2","CE1","CE2","CZ"],
        "TRP":["CG","CD1","CD2","CE2","CE3","CZ2","CZ3","CH2"]
    }
    if str_acide_amine in aromatiques:
        if atome_str in aromatiques[str_acide_amine]:
            return True
    return False


def atom_typing(structure, for_models=None): # Typing et valeurs SASA
    """Affecte la surface accessible et le typing des atomes.

    Argument :
        model : Bio.PDB.Model
        for_models : liste
            Contient les indices des modèles sur lequel
            on veut travailler

    Renvoi : None

    Modifie en place les atomes de notre objet
    """
    if for_models is None:
        for_models = [0]

    # Algorithme pour le Calcul de la surface accessible
    sr = ShrakeRupley()

    # Pour chacun des modèles sur lequel on veut travailler
    for indice, model in enumerate(structure):
        if indice not in for_models:
            continue
        # Pour chacune des chaînes, calcul de la surface accessible des atomes
        for chaine in model:
            # Ajout des valeurs sasa aux atomes
            sr.compute(chaine, level="A")
            # Typing de l'atome, SASA et élément
            for res in chaine:
                res_name = res.get_resname() # Nom du résidu
                for atome in res:
                    a_name = atome.get_name() # Nom de l'atome
                    
                    # Affectation de l'élément
                    if is_Calpha(a_name): # CA
                        atome.element = 'a'
                    elif is_Cbeta(a_name): # CB
                        atome.element = 'b'
                    elif is_CAromatique(res_name,a_name): # C aromatique
                        atome.element = 'A'
                    
                    atome.set_bfactor(atome.sasa) # Affectation du SASA


# ***** Gestion des arguments
# argparse permet de parser les options lors de l'exécution du script
# Nous exigereons d'exécuter notre fichier python en ajoutant en argument
#   obligatoire : le nom du fichier pdb à traiter
#   optionnelle : la localisation du répertoire en sortie (.../PDB_Atyping/)
# python pdb_file [-o] repository_out

args_parser = argparse.ArgumentParser(
    description="""
    Lit un fichier pdb et la renvoi avec les valeurs SASA et le typing de l'élément

    Si l'option [-o] n'est pas spécifié l'emplacement utilisé
        est le même que la où se trouve notre fichier pdb en entrée

    Description :
        La colonne b-factor est remplacée par la valeur SASA
        L'élément change selon ces critères :
        Carbone alpha -> a
        Carbone bétâ -> b
        Carbone aromatique -> A
        On ne modifie pas les éléments des autres atomes.
    """,
    formatter_class=RawTextHelpFormatter
    )

args_parser.add_argument('pdb_file', type=str,help='Chemin vers le fichier pdb') # argument obligatoire
args_parser.add_argument('-o','--output', type=str, # argument optionnel
                        help="Emplacement où sera enregistré notre fichier en sortie")

# Valeurs des arguments passés par l'utilisateur
args = args_parser.parse_args() # Valeurs des arguments passés en argument par l'utilisateur

path_in = args.pdb_file # Nom du chemin vers le fichier pdb
path_out = args.output # Chemin de destination du fichier en sortie

# Vérification de l'existence du fichier pdb
path_exists = os.path.exists(path_in) # renvoi True si le fichier existe
if not path_exists:
    error_message = \
    f"""
    Le fichier pdb dont vous avez spéficié le chemin n'éxiste pas
    path="{path_in}"\n
    """
    sys.exit(error_message)

# PDB en entrée : Emplacement et nom de notre fichier pdb reçu en entrée
path_to_pdb = args.pdb_file.split("/")

pdb_directory = "./" # Correspondra à l'emplacement du répertoire du fichier pdb
pdb_file = "" # nom du fichier pdb complet
pdb_file_name = "" # nom du fichier (sans extension)
pdb_file_extension = "" # extension du fichier

if len(path_to_pdb) > 1:
    # On assigne à pdb_directory le répertoire où se trouve notre fichier
    pdb_directory = "/".join(path_to_pdb[:-1])+"/"

pdb_file_name, pdb_file_extension = os.path.splitext(path_to_pdb[-1]) #nom et extension du fichier pdb
pdb_file = path_to_pdb[-1] # Nom complet (avec extension) de notre fichier

# PDB en sortie : Emplacement et noms de notre répertoire et fichier pdb en sortie
pdb_out_directory_name = "PDB_Atyping" # Nom du répertoire que l'on va créer
pdb_out_directory = os.path.join(pdb_directory, pdb_out_directory_name)

if path_out != None: # | Chemin vers le répertoire à créer
    pdb_out_directory = os.path.join(path_out, pdb_out_directory_name)

pdb_out_file = pdb_file_name + "_Atyping" + pdb_file_extension # nom du fichier pdb à créer
# Chemin complet de notre fichier pdb en sortie
path_out_file = os.path.join(pdb_out_directory, pdb_out_file) # chemin du fichier pdb en sortie


# ***** PDBParser et ShrakeRupley
# Lecture des données PDB à l'aide de PDBParser
parser = PDBParser(QUIET=True) # Lecture Fichier PDB; QUIET=TRUE n'affiche pas de msg d'erreur
structure = parser.get_structure(pdb_file_name, path_in) # PDBParser
work_on_models = [0] # On ne s'intéresse qu'au premier modèle pdb

atom_typing(structure, work_on_models) # SASA et typing
chain_by_models = aux.sections_dict(structure, work_on_models)

model_keys = list(chain_by_models.keys())[0]
sections_to_write_out = aux.str_sections(chain_by_models[model_keys]) + f"{'END':6s}"

# Chemin du répertoire en sortie, le crée s'il n'existe pas
if not os.path.isdir(pdb_out_directory):
    os.mkdir(pdb_out_directory)

# Crée le fichier pdb annoté
with open(path_out_file, "w") as pdb_out:
    pdb_out.write(sections_to_write_out)

# Décommenter pour le typing avec plusieurs modèles
# Crée le fichier pdb annoté
#with open(path_out_file, "w") as pdb_out:
#   sections_to_write_out = ""
#    # Dans le cas ou il n'y a qu'un seul modèle
#    if len(chain_by_models) == 1:
#        model_keys = list(chain_by_models.keys())[0]
#        sections_to_write_out = aux.str_sections(chain_by_models[model_keys])
#
#    else: # Plusieurs modeles # Ajout de la section model
#        for model_keys,model in chain_by_models.items():
#            sections_to_write_out += aux.str_sections(model,model_keys)
#    
#    end_section = f"{'END':6s}"
#    sections_to_write_out += end_section
#
#    # Ecriture dans notre fichier typé.pdb
#    pdb_out.write(sections_to_write_out)

print(f"Fichier '{pdb_out_file}' créé dans : '{pdb_out_directory}'")