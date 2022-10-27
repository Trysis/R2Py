# R2Py

PatchSearch to python

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


positional arguments:
  pdb_file              Chemin vers le fichier pdb

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT


---- Exemple ----

Hiérarchie :
"." Emplacement actuel
|-> "PDB/" répertoire
  |-> "1AWR.pdb" fichier pdb
  |-> "Documents/" répertoire

I.
Dans le bash :
python Atom_Typing.py PDB/1AWR.pdb

Hiérarchie :
"." Emplacement actuel
|-> "PDB/" répertoire
  |-> "1AWR.pdb" fichier pdb
  |-> "Documents/" répertoire
  |-> "PDB_Atyping/" répertoire
    |-> "1AWR_ATyping.pdb" fichier pdb annoté


II. Spécifier le chemin sortant

python Atom_Typing.py -o PDB/Documents PDB/1AWR.pdb

-> Crée le répertoire "./PDB_ATyping/"
-> Crée le fichier 1AWR_ATyping.pdb à "./PDB_ATyping/1AWR_ATyping.pdb"

hiérarchie :
"." Emplacement actuel
|-> "PDB/" répertoire
  |-> "1AWR.pdb" fichier pdb
  |-> "Documents/" répertoire
    |-> "PDB_Atyping/" répertoire
      |-> "1AWR_ATyping.pdb" fichier pdb annoté

## getBindingSite

usage: getBindingSite.py [-h] [-t TCHAIN] [-l LCHAIN] [-c CUTOFF] [-o OUTPUT] [-r] [-a] pdb_file

    Prend en argument un nom de fichier pdb, des chaînes cibles et le ligand
    Et renvoi les atomes (pour le moment les carbones alpha) qui sont à une
    distance inférieur à 5 Angström de notre ligand


positional arguments:
  pdb_file              Chemin vers le fichier pdb

options:
  -h, --help            show this help message and exit
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