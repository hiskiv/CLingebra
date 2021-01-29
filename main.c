#include "./src/CLingebra.h"
#include<stdio.h>

int main(){
    matrix A,B;
    readMatrix(&A);
    displayMatrix(getInverseMatrix(A));
    return 0;
}