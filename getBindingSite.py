# coding=utf-8

import sys
import os # gestion des fichiers
import argparse # gestion des arguments
from argparse import RawTextHelpFormatter
from Bio.PDB import PDBParser # Lecture fichier PDB
from Bio.PDB import NeighborSearch # Algorithme de NeighborSearch

# Import nos fichiers
import auxiliaires as aux

def neighbor_atoms(ligand_atoms, target_chain_atoms, cut=5, res_mode=False, calpha_mode=False):
    """Renvoi la liste des target atomes en contact avec le ligand. 

    Argument :
        ligand_atoms : List(Bio.PDB.Atom.Atom)
            Liste contenant les atomes du ligand
        target_chain_atoms : List(Bio.PDB.Atom.Atom)
            Liste contenalt les atomes d'une chaîne cible
        cut : int
            Valeur du cutoff en Angstrom
        res_mode : bool
            Mode pour spécifier si l'on veut garder tous les atomes du résidu
            dès qu'il y a contact entre l'un des atomes du résidu avec le ligand
        calpha_mode : bool
            Mode pour spécifier si l'on ne veut conserver que les carbones alpha
            dès qu'il y a contact entre l'un des atomes du résidu avec le ligand

    Renvoi : List(Bio.PDB.Atom.Atom)
        Renvoi la liste des atomes en contacts

    """

    neighbor = NeighborSearch(target_chain_atoms)
    close_atoms = None
    if res_mode:  # On conserve tous les atomes en contact avec le ligand
        close_atoms = {
            close_atm 
            for atome in ligand_atoms 
            for close_res in neighbor.search(atome.coord, cutoff, 'R')
            for close_atm in close_res
        }
    elif calpha_mode:  # On conserve tous les carbones alpha dont le résidu est en contact
        close_atoms = {
            close_atm
            for atome in ligand_atoms
            for close_res in neighbor.search(atome.coord, cutoff, 'R')
            for close_atm in close_res
            if close_atm.name == "CA"
        }
    else:  # On conserve uniquement les atomes en contact
        close_atoms = {
            close_atm 
            for atome in ligand_atoms 
            for close_atm in neighbor.search(atome.coord, cutoff, 'A')
        }
    # Tri des atomes proches par leur numéro d'atome
    close_atoms = sorted(close_atoms, key=lambda item: item.get_serial_number())
    return close_atoms


# Gestion des arguments passés par l'utilisateur
args_parser = argparse.ArgumentParser(
    description="""
    Prend en argument un nom de fichier pdb, des chaînes cibles et le ligand
    Et renvoi les atomes (pour le moment les carbones alpha) qui sont à une
    distance inférieur à 5 Angström de notre ligand
    """,
    formatter_class=RawTextHelpFormatter
    )

args_parser.add_argument('pdb_file', type=str,
                        help='Chemin vers le fichier pdb') # argument obligatoire

args_parser.add_argument('-t','--tchain', type=str, # argument optionnel
                        help="Chaînes cibles pour l'évaluation des distances avec le ligand")

args_parser.add_argument('-l','--lchain', type=str, # argument "optionnel"
                        help="Chaîne du ligand dont on évaluera la distance avec nos chaînes cibles")

args_parser.add_argument('-c','--cutoff', type=float, default=5, # argument optionnel,
                        help="Distance en Angström à laquelle on considère"
                        + " des atomes en contactes (compris entre 4 et 12), vaut 5 par défaut")

args_parser.add_argument('-o','--output', type=float, # argument optionnel
                        help="Chemin vers lequel nous enregistrons nos résultats")

args_parser.add_argument('-r','--residu',action='store_true',
                        help="Option pour afficher les atomes des résidus en contact avec le ligand")

args_parser.add_argument('-a','--calpha',action='store_true',
                        help="Option pour afficher les carbones alpha des résidus en contact avec le ligand")

# On récupère les arguments passés par l'utilisateur
args = args_parser.parse_args() # Récupération des arguments

path_in = args.pdb_file # Nom du chemin vers le fichier pdb
path_out = args.output # Chemin de destination du fichier en sortie
tchains = args.tchain # Chaînes cible 
lchains = args.lchain # Ligand cible
cutoff = args.cutoff # Distance à laquelle on considère un contact

residu = args.residu # Option
calpha = args.calpha # Option

# Vérification de l'existence du fichier pdb
path_exists = os.path.exists(path_in) # renvoi True si le fichier existe
if not path_exists:
    error_message = \
    f"""
    Le fichier pdb dont vous avez spéficié le chemin n'éxiste pas :
    path="{path_in}"\n
    """
    sys.exit(error_message)

# Vérifie si l'utilisateur a passé en argument une valeur
# pour la chaîne du ligand
if lchains == None: 
    error_message = \
    f"""
    Veuillez spécifier la chaîne du ligand option [-l].\n
    """
    sys.exit(error_message)

# Cutoff -> Vérifie que le cuftoff est raisonnable
if cutoff < 4 or cutoff > 12:
    error_message = \
    f"""
    Le cutoff='{cutoff}' spécifié n'est pas compris entre 4 et 12.\n
    """
    sys.exit(error_message)

# PDB en entrée : Emplacement et nom de notre fichier pdb reçu en entrée
path_to_pdb = path_in.split("/") # Liste des noms séparés par le slash

pdb_directory = "./" # Correspondra à l'emplacement du répertoire du fichier pdb
pdb_file = "" # nom du fichier pdb complet
pdb_file_name = "" # nom du fichier (sans extension)
pdb_file_extension = "" # extension du fichier
pdb_file = path_to_pdb[-1] # Nom complet (avec extension) de notre fichier
pdb_file_name, pdb_file_extension = os.path.splitext(path_to_pdb[-1]) # nom et extension du fichier pdb

if len(path_to_pdb) > 1: # Assigne à pdb_directory le répertoire où se trouve notre fichier
    pdb_directory = "/".join(path_to_pdb[:-1])+"/"

# PDB en sortie : Emplacement et noms de notre répertoire et fichier pdb en sortie
pdb_out_directory_name = "PDB_Binding" # Nom du répertoire que l'on va créer
pdb_out_directory = os.path.join(pdb_directory,pdb_out_directory_name)

if path_out != None: # | Chemin vers le répertoire à créer
    pdb_out_directory = os.path.join(path_out,pdb_out_directory_name)

### PDBParser
# Lecture des données PDB à l'aide de PDBParser

parser = PDBParser(QUIET=True) # Lecture Fichier PDB; QUIET=TRUE n'affiche pas de msg d'erreur
structure = parser.get_structure(pdb_file_name, path_in) # PDBParser

# Récupération des chaînes passées en argument
pdbchains_set = {chaine_id.id for chaine_id in structure[0]} # chaines du fichier pdb

tchains_set = None # Target chains
if tchains is not None: # Si arguments passés par l'utilisateur
    tchains_set = set(tchains.replace(","," ").strip().split()) 
else: # Sinon correspond à toutes les chaînes de fichier pdb en entrée
    tchains_set = pdbchains_set

lchains_set = set(lchains.strip()) # Ligand chain

# On vérifie que les chaînes passées en arguments existent dans notre fichier pdb
tl_chains_set = tchains_set | lchains_set # Union des chaines de tchains et lchains
diff_chains = tl_chains_set - pdbchains_set # Chaines non existantes dans notre fichier pdb
tchains_set = tchains_set - lchains_set # Pour enlever la même chaîne (ex: A-A)

if len(tchains_set) == 0:
    error_message = \
        f"""
        Vous ne pouvez pas uniquement comparer les même chaînes\n
        """
    sys.exit(error_message)

if len(diff_chains)>0: # Si l'utilisateur spécifie des chaînes non existantes
    error_message = \
        f"""
        Les chaînes {diff_chains} n'existent pas\n
        """
    sys.exit(error_message)

if len(lchains_set)>1:
    error_message = \
        f"""
        Une seule chaîne doit être spécifié pour le ligand.
        """
    sys.exit(error_message)

# Création de notre répertoire pour la conservation du fichier pdb sortant
if not os.path.isdir(pdb_out_directory): # verifie si le répertoire n'existe pas déjà
    os.mkdir(pdb_out_directory)

# ******* getBindingSite
num_modele = 0
latoms = [atome for atome in structure[num_modele][lchains].get_atoms()]  # Atomes du ligand

# On écrira un fichier pour chaque chaînes cibles avec la chaîne du ligand
for tchain in tchains_set:
    # Target atoms de la chaîne tchain
    tatoms_chain = [atome for atome in structure[num_modele][tchain].get_atoms()]
    neighbor_atom_from_tchain = neighbor_atoms(latoms, tatoms_chain, cut=cutoff, res_mode=residu, calpha_mode=calpha)
    
    # pdb_out_file, nom du fichier pdb à créer
    pdb_out_file = pdb_file_name + "_" + tchain + "-" + lchains + pdb_file_extension
    path_out_file = os.path.join(pdb_out_directory, pdb_out_file) # chemin du fichier pdb en sortie
    
    str_section = "".join([aux.atom_section_pdb(atm_objet)+"\n" for atm_objet in neighbor_atom_from_tchain])
    with open(path_out_file, "w") as pdb_out:
        pdb_out.write(str_section)
