#ifndef CLINGEBRA_H_
#define CLINGEBRA_H_

#define MAX_DIM 10

typedef struct _vector{
    int dim;
    double cod[MAX_DIM];
}vector;

typedef struct _matrix{
    int N,M;
    double element[MAX_DIM][MAX_DIM];
}matrix;

//basic vector arithmetic

vector vectorInverse(vector);
vector vectorAddition(vector,vector);
vector vectorSubtraction(vector,vector);
vector vectorScalar(vector,double);
double innerProduct(vector,vector);
vector exteriorProduct(vector,vector);

//basic matrix arithmetic

matrix matrixAddition(matrix,matrix);
matrix matrixMultiplication(matrix,matrix);
matrix matrixTranspose(matrix);
matrix getIdentityMatrix(int);

//matrix and vector discrete operations

matrix vectorsToMatrix(const vector*,int);
vector* matrixToVectors(const matrix);
matrix matrixHorizontallyConcating(matrix,matrix);
matrix matrixVerticallyConcating(matrix,matrix);
matrix matrixHorizontallySlice(matrix,int,int);
matrix matrixVerticallySlice(matrix,int,int);

//basic matrix and vector applications

matrix gaussianElimination(matrix,int);
matrix getInverseMatrix(matrix);
int getRank(matrix);
double getDeterminant(matrix);
double* restrictedGetEigenvalues(matrix);

//handy API for debug and visualizing

void displayVector(vector);
void displayMatrix(matrix);
void readVector(vector*);
void readMatrix(matrix*);

#endif