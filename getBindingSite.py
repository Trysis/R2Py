# coding=utf-8

import sys
import os # gestion des fichiers
import argparse # gestion des arguments
from argparse import RawTextHelpFormatter
from Bio.PDB import PDBParser # Lecture fichier PDB
from Bio.PDB import NeighborSearch # Algorithme de NeighborSearch

# Ecriture des sections dans un fichier pdb
def atom_section_pdb(atom_object):
    """
    Argument:
        objet : atom_object
            Bio.PDB.Atom

    Description:
        Renvoi une ligne de section ATOM (pdb)
        correspondant à l'objet atome passé
        en argument.
    
    Renvoi:
        str : La section (pdb) de l'atome passé en argument
    """
    atom_type = "ATOM"
    # On récupère les arguments qui compose la section ATOM/HETATM
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
    
    if res_id[0] != " ": # Hétéro-atome
        atom_type = "HETATM"

    pdb_atom_str = ""
    pdb_atom_str += f"{atom_type:6s}{serial:5d} {atom_name:^4s}{alt_location:1s}{res_name:3s} {chain_name:1s}"
    pdb_atom_str += f"{res_id[1]:4d}{res_code:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}"
    pdb_atom_str += f"{b_factor:>6.0f}          {element_symbol:>2s}{charge:2s}"

    return pdb_atom_str

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

args_parser.add_argument('-t','--tchain', type=str, # argument optionnel
help="Chaînes cibles pour l'évaluation des distances avec le ligand")

args_parser.add_argument('-l','--lchain', type=str, # argument "optionnel"
help="Chaîne du ligand dont on évaluera la distance avec nos chaînes cibles")

args_parser.add_argument('-c','--cutoff', type=float, # argument optionnel,
default=5, help="Distance en Angström à laquelle on considère"
            +" des atomes en contactes (compris entre 4 et 12), vaut 5 par défaut")

args_parser.add_argument('-o','--output', type=float, # argument optionnel
help="Chemin vers lequel nous enregistrons nos résultats")

args_parser.add_argument('-r','--residu',action='store_true',
help="Option pour afficher les atomes des résidus en contact avec le ligand")

args_parser.add_argument('-a','--calpha',action='store_true',
help="Option pour afficher les carbones alpha des résidus en contact avec le ligand")


# Valeurs des arguments passés par l'utilisateur
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
    Veuillez spécifier la chaîne du ligand option [-l].
    """
    sys.exit(error_message)

if cutoff < 4 or cutoff > 12:
    error_message = \
    f"""
    Le cutoff='{cutoff}' spécifié n'est pas compris entre 4 et 12.
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

### PDBParser
# Lecture des données PDB à l'aide de PDBParser

parser = PDBParser(QUIET=True) # Lecture Fichier PDB; QUIET=TRUE n'affiche pas de msg d'erreur
structure = parser.get_structure(pdb_file_name, path_in) # PDBParser

# Récupération des chaînes passées en argument
pdbchains_set = {chaine_id.id for chaine_id in structure[0]} # chaines du fichier pdb
tchains_set = None
if tchains is not None:
    tchains_set = set(tchains.replace(","," ").split()) # Target chains
else:
    tchains_set = pdbchains_set

lchains_set = set(lchains.strip()) # Ligand chain

# On vérifie dans un premier temps que les chaînes passées en arguments
# sont des chaînes qui existent dans notre fichier pdb
tl_chains_set = tchains_set | lchains_set # Union des chaines de tchains et lchains
diff_chains = tl_chains_set - pdbchains_set # Chaines non existantes dans notre fichier pdb

# Si l'utilisateur a spécifié des chaînes non existantes
if len(diff_chains)>0:
    error_message = \
        f"""
        La chaîne {diff_chains} n'existe pas\n
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

latoms = [atome for atome in structure[0][lchains].get_atoms()]

for chain_name in tchains_set:
    atomes_from_chain = [atome for atome in structure[0][chain_name].get_atoms()]
    neighbor = NeighborSearch(atomes_from_chain)
    close_atoms = None
    
    if residu :
        close_atoms = sorted(
        {close_atm 
        for atome in latoms 
        for close_res in neighbor.search(atome.coord, cutoff, 'R')
        for close_atm in close_res}
        ,
        key=lambda item: item.get_serial_number()
        )
    elif calpha:
        close_atoms = sorted(
        {close_atm
        for atome in latoms
        for close_res in neighbor.search(atome.coord, cutoff, 'R')
        for close_atm in close_res
        if close_atm.name == "CA"},
        key=lambda item: item.get_serial_number()
        )
    else: #Atomes en contact
        close_atoms = sorted(
            {close_atm for atome in latoms for close_atm in neighbor.search(atome.coord, cutoff, 'A')},
            key=lambda item: item.get_serial_number()
            )
    str_section = "".join([atom_section_pdb(atm_objet)+"\n" for atm_objet in close_atoms])
    
    # pdb_out_file, nom du fichier pdb à créer
    pdb_out_file = pdb_file_name + "_" + chain_name + ":" + lchains + pdb_file_extension
    # Chemin complet
    path_out_file = os.path.join(pdb_out_directory, pdb_out_file) # chemin du fichier pdb en sortie

    with open(path_out_file,"w") as pdb_out:
        pdb_out.write(str_section)
