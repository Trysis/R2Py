# coding=utf-8

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

def ter_section_pdb(chain_id,res,atome):
    """
    Arguments:
        str : chain_id
            Nom de la chaine de l'atome
        objet : res
            Bio.PDB.Residue (residue de l'atome)
        objet : atome
            Bio.PDB.Atom (atome)

    Description:
        Renvoi une ligne de section TER en
        prenant en compte le nom de la chaîne
        de l'atome passé en argument et son résidue
    
    Renvoi:
        str : Renvoi une ligne de la section TER
    """
    atome_serial = atome.get_serial_number()
    ter_section = f"{'TER':6s}{atome_serial+1:>5d}      "
    ter_section += f"{res.get_resname():3s} "
    ter_section += f"{chain_id:1s}{res.id[1]:4d}{res.id[0]:1s}\n"

    return ter_section


def str_sections(model_section, model_key=-1):
    """
    Arguments:
        dictionnaire : model_section
            Dictionnaire qui contient pour chaque
            clés 'chaîne' les sections atomes de la chaîne (str)
        int : model_key
            numéro du modèle dans le cas où l'on veut
            ajouter les sections MODEL/ENDMDL

    Description:
        Renvoi les sections atomes et/ou hetatm d'un modèle
    
    Renvoi:
        str : Les sections atomes/hetatm d'un modele
    """
    # Dans le cas ou il n'y a qu'un seul modèle
    if model_key == -1:
        return "".join(model_section.values())

    else: # Ajout de la section model correspondant
        model_section = f"{'MODEL':6s}    {model_key:>4d}\n"
        atm_hetatm_section = "".join(model_section.values())
        endmodel_section = f"{'ENDMDL':6s}\n"
        return model_section + atm_hetatm_section + endmodel_section
