# coding=utf-8

# Ecriture des sections dans un fichier pdb
def atom_section_pdb(atom_object):
    """Renvoi une ligne de section ATOM ou HETATM qui respecte le format pdb.
        La ligne renvoyée est issu des informations de l'objet atome passé en argument.

    Arguments :
        atom_object : Bio.PDB.Atom.Atom
    
    Renvoi : str
        La section pdb ATOM/HETATM de lobjet passé en argument
    """

    # Arguments qui compose la section ATOM/HETATM
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
    
    if res_id[0] != " ": # Hétéro-atome
        atom_type = "HETATM"

    # Écriture sous le bon format dans une chaîne de caractère
    pdb_atom_str = ""
    pdb_atom_str += f"{atom_type:6s}{serial:5d} {atom_name:^4s}{alt_location:1s}{res_name:3s} {chain_name:1s}"
    pdb_atom_str += f"{res_id[1]:4d}{res_code:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}"
    pdb_atom_str += f"{b_factor:>6.0f}          {element_symbol:>2s}{charge:2s}"

    return pdb_atom_str


def ter_section_pdb(chain_id, res, atome):
    """Renvoi une ligne de section TER qui respecte le format pdb.
        La ligne renvoyée prend en compte le nom de la chaîne
        de l'atome passé en argument et le résidu associé.

    Arguments :
        chain_id : str
            Nom de la chaine de l'atome
        res : objet -> Bio.PDB.Residue.Residue
        atome : objet -> Bio.PDB.Atom.Atom

    Renvoi : str
        Renvoi une ligne de la section TER
    """

    atome_serial = atome.get_serial_number()
    ter_section = f"{'TER':6s}{atome_serial+1:>5d}      "
    ter_section += f"{res.get_resname():3s} "
    ter_section += f"{chain_id:1s}{res.id[1]:4d}{res.id[0]:1s}\n"

    return ter_section


### Ecriture des sections à mettre dans le fichier pdb
# Dictionnaire de modèles
# Ecrit au sein de chaque modèles les sections à écrire dans le fichier pdb
#   chain_by_models[model][chaine] = sections atome/hetatm de la chaîne 'chaine'


def sections_dict(structure, for_models=None):
    """Renvoi un dictionnaire contenant pour chaque modèle ses sections pdb à afficher.

    Arguments :
        structure : objet -> Bio.PDB.Structure.Structure
        for_models : liste
            Contient les indices des modèles sur lequel
            on veut travailler

    Renvoi : str
        Renvoi une ligne de la section TER
    """

    if for_models is None:
        for_models = [0]

    chain_by_models = {}
    for indice, model in enumerate(structure):
        if indice not in for_models:
            continue
        last_atm = None
        last_res = None
        hetatm_list = ""

        chain_by_models[model.serial_num] = dict()
        for chaine in model:
            atoms_section_in_chain = ""
            for atome in chaine.get_atoms():
                res = atome.get_parent()
                if res.id[0] == " ": # Atome
                    atoms_section_in_chain += f"{atom_section_pdb(atome)}\n"
                    last_atm = atome
                    last_res = atome.get_parent()
                else: # Hétéro-Atome
                    hetatm_list += f"{atom_section_pdb(atome)}\n"
            # Lorsque l'on passe d'une chaîne à l'autre on écrit la section TER
            chain_by_models[model.serial_num][chaine.id] = \
                atoms_section_in_chain + ter_section_pdb(chaine.id, last_res, last_atm)

        chain_by_models[model.serial_num]["HETATM"] = hetatm_list
    return chain_by_models


def str_sections(model_section, model_key=-1):
    """Renvoi les sections atomes et/ou hetatm d'un modèle.

    Arguments :
        model_section : dictionnaire
            Dictionnaire qui contient pour chaque
            clés 'chaîne' les sections atomes de la chaîne (str)
        model_key : int
            numéro du modèle dans le cas où l'on veut
            ajouter les sections MODEL/ENDMDL        

    Renvoi : str
        Les sections atomes/hetatm d'un modele
    """

    # Dans le cas ou il n'y a qu'un seul modèle
    if model_key == -1:
        return "".join(model_section.values())

    else: # Ajout de la section model correspondant
        model_section = f"{'MODEL':6s}    {model_key:>4d}\n"
        atm_hetatm_section = "".join(model_section.values())
        endmodel_section = f"{'ENDMDL':6s}\n"
        return model_section + atm_hetatm_section + endmodel_section
