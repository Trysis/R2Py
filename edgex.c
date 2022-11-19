#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

//gcc -shared edgex.c -o edgex.so

double dist3D(double x1, double y1, double z1, double x2, double y2, double z2){
    double dx = (x1 - x2);
    double dy = (y1 - y2);
    double dz = (z1 - z2);
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void vertex(uint32_t ** vertex, char * X_element, char * Y_element, size_t XYsize[2], size_t * vertex_size){
    size_t compteur = 0;
    for(size_t i = 0; i < XYsize[0]; i++){
        for(size_t j = 0; j < XYsize[1]; j++){
            if(X_element[i] == Y_element[j]){
                vertex[0][compteur] = i;
                vertex[1][compteur] = j;
                compteur += 1;
            }
        }
    }
    *vertex_size = compteur;
}

void edge(uint32_t ** edge, uint32_t ** indexXY, double ** coordX, double ** coordY, size_t nrow_N_M[3], size_t * edge_size){
    size_t nrow = nrow_N_M[0];
    size_t N = nrow_N_M[1];
    size_t M = nrow_N_M[2];
    double seuil = 1;

    size_t compteur = 0;
    for(size_t i=0; i<nrow; i++){
        for(size_t j=i+1; j<nrow; j++){
            if(indexXY[0][i] == indexXY[0][j])continue;
            if(indexXY[1][i] == indexXY[1][j])continue;

            // coordonnee dans X
            double x1 = coordX[0][i];
            double y1 = coordX[1][i];
            double z1 = coordX[2][i];

            double x1_prm = coordX[0][j];
            double y1_prm = coordX[1][j];
            double z1_prm = coordX[2][j];

            // coordonnee dans Y
            double x2 = coordY[0][i];
            double y2 = coordY[1][i];
            double z2 = coordY[2][i];

            double x2_prm = coordY[0][j];
            double y2_prm = coordY[1][j];
            double z2_prm = coordY[2][j];

            //Calcul des distances
            double d1 = dist3D(x1, y1, z1, x1_prm, y1_prm, z1_prm);
            double d2 = dist3D(x2, y2, z2, x2_prm, y2_prm, z2_prm);
            
            if(abs(d1 - d2) <= seuil){
                edge[0][compteur] = i;
                edge[1][compteur] = j;
                compteur += 1;
            }
        }
    }
    *edge_size = compteur;
}

int ** edge2(int * index_X, int * index_Y, double ** coordX, double ** coordY, int nrow_coord, int N, int M, size_t * nrow_edge){
    int ** edge_matrix = malloc(sizeof(int *) * 2); // edge matrix
    assert(edge_matrix != NULL);

    size_t size = (size_t)(nrow_coord * nrow_coord);
    
    //Allocation mémoire de chaques colonnes
    for(int i=0; i<2; i++){
        edge_matrix[i] = malloc(sizeof(int) * size);
        assert(edge_matrix[i] != NULL);
    }

    int seuil = 1; // Seuil en Angstrom
    size_t nrow_e = 0;
    for(int i=0; i<nrow_coord; i++){
        for(int j=i+1; j<nrow_coord; j++){
            if(*nrow_edge == (size_t) 48043500)printf("j=%d\n", j);
            if(index_X[i] == index_X[j])continue;
            if(index_Y[i] == index_Y[j])continue;

            // coordonnee dans X
            double x1 = coordX[0][i];
            double y1 = coordX[1][i];
            double z1 = coordX[2][i];

            double x1_prm = coordX[0][j];
            double y1_prm = coordX[1][j];
            double z1_prm = coordX[2][j];

            // coordonnee dans Y
            double x2 = coordY[0][i];
            double y2 = coordY[1][i];
            double z2 = coordY[2][i];

            double x2_prm = coordY[0][j];
            double y2_prm = coordY[1][j];
            double z2_prm = coordY[2][j];

            //Calcul des distances
            double d1 = dist3D(x1, y1, z1, x1_prm, y1_prm, z1_prm);
            double d2 = dist3D(x2, y2, z2, x2_prm, y2_prm, z2_prm);
            
            if(abs(d1 - d2) <= seuil){
                /*
                int md_idx = index_X[i]*N + index_Y[i];
                int md_idx_prm = index_X[j]*N + index_Y[j];
                
                printf("N = %d; M = %d\n", N, M);
                printf("i:%d [%f, %f, %f]; j:%d [%f, %f, %f]\n",
                    index_X[i], x1, y1, z1,
                    index_Y[i], x2, y2, z2);
                printf("i':%d [%f, %f, %f]; j':%d [%f, %f, %f]\n",
                    index_X[j], x1_prm, y1_prm, z1_prm,
                    index_Y[j], x2_prm, y2_prm, z2_prm);
                
                printf("md_idx = %d; md_idx_prm = %d\n", md_idx, md_idx_prm);
                printf("d1 = %f; d2 = %f\n\n", d1, d2);
                */

                edge_matrix[0][*nrow_edge] = i;
                edge_matrix[1][*nrow_edge] = j;

                nrow_e += (size_t) 1;
                if(nrow_e >= size){
                    size += (size_t)(nrow_coord * nrow_coord);
                    for(int i=0; i<2; i++){
                            void * v = realloc(edge_matrix[i], sizeof(int) * size);
                            assert(v != NULL);
                    }
                }
            }
        }
    }
    *nrow_edge = nrow_e;
    printf("hey");
    // Reallocation de la memoire au nombre de lignes attribuee
    for(int i=0; i<2; i++){
        void * v = (int *) realloc(edge_matrix[i], sizeof(int) * (size_t)(*nrow_edge));
        assert(v != NULL);
    }
    return edge_matrix;
}
