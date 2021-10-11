#include "./src/CLingebra.h"
#include<stdio.h>

int main(){
    matrix A,B;
    readMatrix(&A);
    double *eigenvalues = restrictedGetEigenvalues(A);
    for(int i=1;i<=eigenvalues[0];i++){
        printf("%.8lf\n", eigenvalues[i]);
    }
    return 0;
}