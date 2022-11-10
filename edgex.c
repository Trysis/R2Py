#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//gcc -shared edgex.c -o edgex.so

void edge(int ** elem, int x){
    // x nombre de lignes
    // y nombre de colonnes
    for(int i = 0; i<x; i++){
        printf("x=%d y=%d\n", elem[0][i], elem[1][i]);
    }
}

int dist3D(double x1, double y1, double z1, double x2, double y2, double z2){
    double dx = (x1 - x2);
    double dy = (y1 - y2);
    double dz = (z1 - z2);
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void edge2(int * elem_idx1, int * elem_idx2, double ** coord1, double ** coord2, int nrow){
    int ** edge_matrix = malloc(sizeof(int*) * 2);
    int mat_row = nrow * nrow;
    printf("hey\n")
    //Allocation mémoire des chaques colonnes
    for(int i=0; i<2; i++)edge_matrix[i] = malloc(sizeof(int) * mat_row);
    for(int i=0; i<nrow; i++){
        for(int j=i+1; j<nrow; j++){
            if(elem_idx1[i] == elem_idx1[j])continue;
            if(elem_idx2[i] == elem_idx2[j])continue;

            double x1 = coord1[i][0];
            double y1 = coord1[i][1];
            double z1 = coord1[i][2];

            double x2 = coord2[i][0];
            double y2 = coord2[i][1];
            double z2 = coord2[i][2];

            double x1_prm = coord1[j][0];
            double y1_prm = coord1[j][1];
            double z1_prm = coord1[j][2];

            double x2_prm = coord2[j][0];
            double y2_prm = coord2[j][1];
            double z2_prm = coord2[j][2];

            double d1 = dist3D(x1, y1, z1, x1_prm, y1_prm, z1_prm);
            double d2 = dist3D(x2, y2, z2, x2_prm, y2_prm, z2_prm);
            printf("[%d,%d]:%f ; [%d,%d]:%f\n", i, j, d1, i, j, d2);
            return;
        }
    }
}
