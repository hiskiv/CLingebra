#include "./src/CLingebra.h"
#include<stdio.h>

int main(){
    matrix A,B;
    readMatrix(&A),readMatrix(&B);

    displayMatrix(matrixHorizontallyConcating(A,B));
    return 0;
}