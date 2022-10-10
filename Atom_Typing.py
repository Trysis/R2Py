# coding=utf-8

import sys
import os # gestion des fichiers
import argparse # gestion des arguments
from argparse import RawTextHelpFormatter

from Bio.PDB import PDBParser # Lecture fichier PDB
from Bio.PDB.SASA import ShrakeRupley # Algorithme de calcul de la surface accessible

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

# Fonctions

def is_Calpha(str_carbone):
    """Renvoi True si str_carbone est un carbone alpha"""
    return "CA" == str_carbone


def is_Cbeta(str_carbone):
    """Renvoi True sir str_carbone est un carbone bétâ"""
    return "CB" == str_carbone


def is_Carbone_R(str_carbone):
    """Renvoi True si str_carbone est un carbone de la chaîne latérale"""
    if len(str_carbone) <= 1: # Verifie que nous ne sommes pas sur le carbone COO
        return False
    return "C" == str_carbone[0]


def is_aromatique(str_acide_amine):
    """Renvoi True sir str_acide_amine est un résidu aromatique"""
    aromatiques_valides = ["HIS","PHE","TYR","TRP"]
    return str_acide_amine in aromatiques_valides

def atom_to_pdb(atom_object):
    atom_type = "ATOM"
    serial = atom_object.get_serial_number()
    atom_name = atom_object.get_name()
    alt_location = atom_object.get_altloc()
    res_name = atom_object.get_parent().get_resname()
    _,_,chain_name,res_id,_ = atom_object.get_full_id()
    
    res_code = res_id[2]
    x, y, z = atom_object.get_coord()
    occupancy = atom_object.get_occupancy()
    b_factor = atom_object.get_bfactor()
    element_symbol = atom_object.element
    charge = " "
    
    if res_id[0] != " ":
        atom_type = "HETATM"

    pdb_atom_str = ""
    pdb_atom_str += f"{atom_type:6s}{serial:5d} {atom_name:^4s}{alt_location:1s}{res_name:3s} {chain_name:1s}"
    pdb_atom_str += f"{res_id[1]:4d}{res_code:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}"
    pdb_atom_str += f"{b_factor:>6.0f}          {element_symbol:>2s}{charge:2s}"

    return pdb_atom_str

def ter_section_pdb(chain_id,res,atome):
    atome_serial = atome.get_serial_number()
    ter_section = f"{'TER':6s}{atome_serial+1:>5d}      "
    ter_section += f"{res.get_resname():3s} "
    ter_section += f"{chain_id:1s}{res.id[1]:4d}{res.id[0]:1s}\n"

    return ter_section


# Gestion des arguments
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
args_parser.add_argument('-o','--path_out', type=str, # argument optionnel
help="Emplacement où sera enregistré notre fichier en sortie")

args = args_parser.parse_args() # Valeurs des arguments passés en argument par l'utilisateur

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
if args.path_out != None:
    pdb_out_directory = os.path.join(args.path_out,"PDB_Atyping")

pdb_out_file = pdb_file_name + "_Atyping" + pdb_file_extension # nom du fichier pdb à créer

# Chemin complets pour nos fichiers
path_file = args.pdb_file # chemin du fichier pdb en entrée
path_out_file = os.path.join(pdb_out_directory, pdb_out_file) # chemin du fichier pdb en sortie

# Création et vérification fichier/répertoire

# Vérification de l'existence des fichier
path_exists = os.path.exists(path_file) # verifie si le fichier pdb existe
if not path_exists:
    error_message = \
    f"""
    Le fichier pdb dont vous avez spéficié le chemin n'éxiste pas
    path="{path_file}"\n
    """
    sys.exit(error_message)

# Création de notre répertoire pour la conservation du fichier pdb sortant
if not os.path.isdir(pdb_out_directory): # verifie si le répertoire n'existe pas déjà
    os.mkdir(pdb_out_directory)

# Biopython
# PDBParser
# Lecture des données PDB à l'aide de PDBParser

parser = PDBParser(QUIET=True) # Lecture Fichier PDB; QUIET=TRUE n'affiche pas de msg d'erreur
sr = ShrakeRupley() # Calcul de la surface accessible

structure = parser.get_structure(pdb_file_name, path_file) # PDBParser

# Calcul de la surface accessible de chacun des atomes pour chaques chaînes de notre protéine
# Et typing de l'élément
for model in structure:
    # Calcul de la surface accessible
    for chaine in model:
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
                elif is_aromatique(res_name) and is_Carbone_R(a_name): # C aromatique
                    atome.element = 'A'
                
                atome.set_bfactor(atome.sasa) # Affectation du SASA

# Ecriture des sections atomes dans un dictionnaire 
#   de chaîne contenu dans un dictionnaire de modeles
# Soit :
#   chains_from_models[model][chaine] = atomes sections de la chaîne
#
chains_from_models = {}
hetatm_list = ""

for model in structure:
    last_atm = None
    last_res = None

    chains_from_models[model.serial_num] = dict()
    for chaine in model:
        for atome in chaine.get_atoms():
            res = atome.get_parent()
            if res.id[0] == " ": # Atome
                chains_from_models[model.serial_num][chaine.id] = \
                    chains_from_models[model.serial_num].get(chaine.id,"")\
                        + f"{atom_to_pdb(atome)}\n"
                last_atm = atome
                last_res = atome.get_parent()
            else: # Hétéro-Atome
                hetatm_list += f"{atom_to_pdb(atome)}\n"
        # Lorsque l'on passe d'une chaîne à l'autre on écrit la section TER
        chains_from_models[model.serial_num][chaine.id] += \
            ter_section_pdb(chaine.id, last_res, last_atm)


#Création de notre fichier pdb typé
with open(path_out_file, "w") as pdb_out:
    section_to_write_out = ""
    # Dans le cas o il n'y a qu'un seul modèle
    if len(chains_from_models) == 1:
        atomes_section = "".join(chains_from_models[0].values())
        hetatm_section = hetatm_list
        section_to_write_out+= atomes_section+hetatm_section

    else: # Plusieurs modeles # Ajout de la section model
        for key_model in chains_from_models:
            model_section = f"{'MODEL':6s}    {key_model:>4d}\n"
            atomes_section = "".join(chains_from_models[key_model].values())
            hetatm_section = hetatm_list
            endmodel_section = f"{'ENDMDL':6s}\n"
            section_to_write_out += model_section+atomes_section \
                                    + hetatm_section+endmodel_section
    
    end_section = f"{'END':6s}"
    section_to_write_out += end_section
    pdb_out.write(section_to_write_out)

print(f"Fichier créé dans : {pdb_out_directory}")