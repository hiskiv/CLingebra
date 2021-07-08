#include "./src/CLingebra.h"
#include<stdio.h>

int main(){
    matrix A,B;
    readMatrix(&A);
    printf("%lf\n", getDeterminant(A));
    return 0;
}