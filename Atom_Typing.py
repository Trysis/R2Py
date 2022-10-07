import os
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB.SASA import ShrakeRupley

# Variables globales

parser = PDBParser(QUIET=True)# Lecture Fichier PDB; QUIET=TRUE n'affiche pas msg d'erreur
io = PDBIO() # Sauvegarde fichier pdb seulements model/atomes
sr = ShrakeRupley() # Calcul de la surface accessible

atome_pos_record = \
    {
    "ATOM":(1,4),"serial":(7,11),"atome_name":(13,16),
    "a_location":(17,17),"res_name":(18,20),"chain":(22,22),
    "res_num":(23,26),"res_code":(27,27),
    "x":(31,38),"y":(39,46),"z":(47,54),
    "occupancy":(55,60),"temperature":(61,66),
    "segment":(73,76),"element":(77,78)
    } # Position de chacune des colonnes dans un fichier pdb

a_pos_record = {key:slice(val[0]-1,val[1]) for key,val in atome_pos_record.items()}


# Fonctions

def is_Carbone_R(str_carbone): # Verifie si Carbone d'une chaîne latérale
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

# Classe a utiliser pour lors de la creation de notre fichier à l'aide du module
#   PDBIO, permet de spécifier quelles atomes/résidus/chaines/modeles nous voulons garder
#   ici je l'utilise pour qu'il affecte la valeur SASA au b_factor juste avant la creation du fichier
class CASelect(Select):
    def accept_atom(self,atom):
        atom.set_bfactor(atom.sasa)
        return 1
# Main

pdb_directory = "PDB" # repertoire ou se situe le fichier .pdb

pdb_file_name = "2ptn" # nom du fichier
pdb_file_extension = ".pdb" # extension du fichier
pdb_file = pdb_file_name + pdb_file_extension # nom complet du fichier pdb

pdb_output_directory = pdb_directory+"_patch"
pdb_output_file = pdb_file_name + "_patch" + pdb_file_extension # nom du fichier pdb modifié

path_file = os.path.join(pdb_directory, pdb_file) # chaîne "pdb_dicretory/pdb_file"
path_output_file = os.path.join(pdb_output_directory, pdb_output_file)

path_exists = os.path.exists(path_file) # verifie si le fichier pdb existe

print(f"Fichier existe={str(path_exists)}\n > path={path_file} \n")

#BioPDB lecture du fichier pdb
structure = parser.get_structure(pdb_file_name, path_file) # PDBParser
io.set_structure(structure) # nous permettra de créer un fichier pdb

# Calcul de la surface accessible de chacun des atomes pour chaques chaînes de notre protéine
for model in structure:
    for chaine in model:
        sr.compute(chaine, level="A")

# Creation de notre nouveau fichier pdb avec valeur SASA
if not os.path.isdir(pdb_output_directory): # verifie si le fichier pdb existe
    os.mkdir(pdb_output_directory)

io.save(path_output_file,CASelect())
path_exists = os.path.exists(path_output_file) # verifie le chemin du fichier pdb modifie

#Ouverture du fichier pdb en lecture et ecriture pour ajouter les element ('a','b','A')
with open(path_output_file, "r+") as pdb_out:
    offset = 0 # position dans le cadre de lecture
    numero_ligne = 0 # numero de la ligne dans le fichier pdb

    # taille de la section element
    element_length = a_pos_record["element"].stop-a_pos_record["element"].start 
    
    for lines in pdb_out:
        if "ATOM" not in lines[a_pos_record["ATOM"]]: 
            pass

        atome_name = lines[a_pos_record["atome_name"]]# nom de l'atome a notre ligne
        res_name = lines[a_pos_record["res_name"]] # nom du residu a notre ligne

        if is_Calpha(atome_name):
            # On modifie le nom de l'élement à "a" si carbone alpha
            str_to_write = f"{'a':<{element_length}}"
            pdb_out.seek(offset+a_pos_record["element"].start+1)
            pdb_out.write(str_to_write)

        elif is_Cbeta(atome_name):
            # On modifie le nom de l'élement à "b" si carbone beta
            str_to_write = f"{'b':<{element_length}}"
            pdb_out.seek(offset+a_pos_record["element"].start+1)
            pdb_out.write(str_to_write)

        elif is_aromatique(res_name):
            if is_Carbone_R(atome_name):
                # Si residu aromatique, et carbone R, element devient "A"
                str_to_write = f"{'A':<{element_length}}"
                pdb_out.seek(offset+a_pos_record["element"].start+1)
                pdb_out.write(str_to_write)
        
        offset += len(lines) # Position de notre cadre de lecture
        pdb_out.seek(offset) # Necessaire si l'on veut maintenir la position de notre ligne

print("job done")