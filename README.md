# R2Py

PatchSearch en python

## Atom_Typing.py

Ce script python permet l'annotation des valeurs SASA de chacuns des atomes, et l'annotation des éléments.

usage: Atom_Typing.py [-h] [-o OUTPUT] pdb_file

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


arguments obligatoires:
  pdb_file              Chemin vers le fichier pdb

options:
  -h, --help            affiche le message d'aide
  -o OUTPUT, --output OUTPUT


---- Exemple ----

I. Utilisation simple

python Atom_Typing.py PDB/1AWR.pdb

-> Crée le fichier "1AWR_ATyping.pdb" à "./PDB/PDB_ATyping/1AWR_ATyping.pdb"

II. Spécifier le chemin sortant

python Atom_Typing.py -o . PDB/1AWR.pdb

-> Crée le répertoire "./PDB_ATyping/"
-> Crée le fichier "1AWR_ATyping.pdb" à "./PDB_ATyping/1AWR_ATyping.pdb"


## getBindingSite

usage: getBindingSite.py [-h] [-t TCHAIN] [-l LCHAIN] [-c CUTOFF] [-o OUTPUT] [-r] [-a] pdb_file

    Prend en argument un nom de fichier pdb, des chaînes cibles et le ligand
    Et renvoi les atomes qui sont à une
    distance inférieur à 5 Angström de notre ligand


arguments obligatoires:
  pdb_file              Chemin vers le fichier pdb

options:
  -h, --help            message d'aide
  -t TCHAIN, --tchain TCHAIN
                        Chaînes cibles pour l'évaluation des distances avec le ligand
  -l LCHAIN, --lchain LCHAIN
                        Chaîne du ligand dont on évaluera la distance avec nos chaînes cibles
  -c CUTOFF, --cutoff CUTOFF
                        Distance en Angström à laquelle on considère des atomes en contactes (compris entre 4 et 12), vaut 5 par défaut
  -o OUTPUT, --output OUTPUT
                        Chemin vers lequel nous enregistrons nos résultats
  -r, --residu          Option pour afficher les atomes des résidus en contact avec le ligand
  -a, --calpha          Option pour afficher les carbones alpha des résidus en contact avec le ligand

---- Exemple ----

I. Atomes en contacts avec les atomes de la chaîne A

python getBindingSite.py -l A PDB/1AWR.pdb

-> Crée un fichier "./PDB/PDB_Binding/1AWR_P-A" (1AWR possède 2 chaînes, A et P)

II. Carbones alpha quand le résidus est en contact avec les atomes de la chaîne A
Attention ! (Écrase le fichier de même nom)

python getBindingSite.py -l A PDB/1AWR.pdb --calpha

-> Crée un fichier "./PDB/PDB_Binding/1AWR_P-A"

III. Affiche tous les atomes d'un résidus lorsqu'il est en contact avec un atome de la chaîne A

python getBindingSite.py -l A PDB/1AWR.pdb --residu

-> Crée un fichier "./PDB/PDB_Binding/1AWR_P-A"
