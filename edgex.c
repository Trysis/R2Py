#include <stdio.h>

void edge(int ** elem, int x){
    // x nombre de lignes
    // y nombre de colonnes
    for(int i = 0; i<x; i++){
        printf("x=%d y=%d\n", elem[0][0], elem[1][1]);
    }
}