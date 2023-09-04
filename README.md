# R2Py  

PatchSearch[1] en python  

1. Rey J, Rasolohery I, Tufféry P, Guyon F, Moroy G. PatchSearch: a web server for off-target protein identification. Nucleic Acids Res. 2019 Jul 2;47(W1):W365-W372. doi: 10.1093/nar/gkz478. PMID: 31131411; PMCID: PMC6602448.


## Packages à installer  

```conda install -c conda-forge biopandas numpy python-igraph pandas```



## Atom_Typing.py  
-----------------  
Ce script python permet l'annotation des valeurs SASA de chacuns des atomes, et l'annotation des éléments.

usage: `python Atom_Typing.py [-h] [-o OUTPUT] pdb_file`  
  
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

__I. Utilisation simple__

```
python Atom_Typing.py PDB/1AWR.pdb
```  

-> Crée le fichier "1AWR_ATyping.pdb" à "./PDB/PDB_ATyping/1AWR_ATyping.pdb"

__II. Spécifier le chemin sortant__

```
python Atom_Typing.py -o . PDB/1AWR.pdb
```  

-> Crée le répertoire "./PDB_ATyping/"
-> Crée le fichier "1AWR_ATyping.pdb" au sein de l'emplacement "./PDB_ATyping/1AWR_ATyping.pdb"



## getBindingSite  
----------------- 
usage: `getBindingSite.py [-h] [-t TCHAIN] [-l LCHAIN] [-c CUTOFF] [-o OUTPUT] [-r] [-a] pdb_file`

    Prend en argument un nom de fichier pdb, des chaînes cibles et le ligand
    Et renvoi les atomes qui sont à une
    distance inférieur à 5 Angström de notre ligand


arguments obligatoires:  
  pdb_file              Chemin vers le fichier pdb  

options:  
  -h, --help &emsp;&emsp;&emsp; message d'aide  
  -t TCHAIN, --tchain TCHAIN  
  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Chaînes cibles pour l'évaluation des distances avec le ligand  
  -l LCHAIN, --lchain LCHAIN  
  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Chaîne du ligand dont on évaluera la distance avec nos chaînes cibles  
  -c CUTOFF, --cutoff CUTOFF  
  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Distance en Angström à laquelle on considère des atomes en contactes (compris entre 4 et 12), vaut 5 par défaut  
  -o OUTPUT, --output OUTPUT  
  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Chemin vers lequel nous enregistrons nos résultats  
  -r, --residu &emsp;&emsp;&emsp; Option pour afficher les atomes des résidus en contact avec le ligand  
  -a, --calpha &emsp;&emsp;&emsp; Option pour afficher les carbones alpha des résidus en contact avec le ligand  

---- Exemple ----  

__I. Atomes en contacts avec les atomes de la chaîne A__

```
python getBindingSite.py -l A PDB/1AWR.pdb
```  

-> Crée un fichier "./PDB/PDB_Binding/1AWR_P-A" (1AWR possède 2 chaînes, A et P)

__II. Carbones alpha quand le résidus est en contact avec les atomes de la chaîne A Attention !__ (Écrase le fichier de même nom)

```
python getBindingSite.py -l A PDB/1AWR.pdb --calpha
```  

-> Crée un fichier "./PDB/PDB_Binding/1AWR_P-A"

__III. Affiche tous les atomes d'un résidus lorsqu'il est en contact avec un atome de la chaîne A__

```
python getBindingSite.py -l A PDB/1AWR.pdb --residu
```  

-> Crée un fichier "./PDB/PDB_Binding/1AWR_P-A"  

## Motifs (Incomplet, il manque le calcul du RMSD)  
-----------------  

compiler vertex_edge.c :
``` 
gcc -shared vertex_edge.c -o vertex_edge.so  
``` 

usage: Motifs.py ligand.pdb target.pdb

-----------------
