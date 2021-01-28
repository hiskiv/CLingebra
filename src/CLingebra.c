#include "CLingebra.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define iszero(a) (fabs(a)<1e-4)

//Get the inverse element of a vector
vector vectorInverse(vector dst){
    for(int i=0;i<dst.dim;i++){
        dst.cod[i]*=-1;
    }
    return dst;
}

//Add two vectors: dst = dst + src
vector vectorAddition(vector dst,vector src){
    if(dst.dim!=src.dim){
        printf("Dimension mismatching!\n");
        return dst;
    }
    for(int i=0;i<dst.dim;i++){
        dst.cod[i]+=src.cod[i];
    }
    return dst;
}

//Minus two vectors: dst = dst - src
vector vectorMinus(vector dst,vector src){
    vector t=vectorInverse(src);
    return vectorAddition(dst,t);
}

//Scalar multiplication of integer k and vector src
vector vectorScalar(vector src,double k){
    for(int i=0;i<src.dim;i++){
        src.cod[i]*=k;
    }
    return src;
}

//Calculate the inner product of vectors src and dst
double innerProduct(vector src,vector dst){
    double ret=0;
    if(src.dim!=dst.dim){
        printf("Dimension mismatching!\n");
        return -1;
    }
    for(int i=0;i<src.dim;i++){
        ret+=src.cod[i]*dst.cod[i];
    }
    return ret;
}

//Calculate the exterior product of vectors src and dst
vector exteriorProduct(vector src,vector dst){}

//Addition of two matrices: dst = dst + src
matrix matrixAddition(matrix dst,matrix src){
    if(dst.N!=src.N||dst.M!=src.M){
        printf("Dimension mismatching!\n");
        return dst;
    }
    for(int i=0;i<dst.N;i++){
        for(int j=0;j<dst.M;j++){
            dst.element[i][j]+=src.element[i][j];
        }
    }
    return dst;
}

//Multiplication of two matrices: dst = dst * src
matrix matrixMultiplication(matrix dst,matrix src){
    if(dst.M!=src.M){
        printf("Multiplication arguments mismatching!\n");
        return dst;
    }
    matrix ret;
    ret.N=dst.N,ret.M=src.M;
    for(int i=0;i<dst.N;i++){
        for(int j=0;j<src.M;j++){
            ret.element[i][j]=0;
            for(int k=0;k<dst.M;k++){
                ret.element[i][j]+=dst.element[i][k]*src.element[k][j];
            }
        }
    }
    return ret;
}

//Transpose the matrix src
matrix matrixTranspose(matrix src){
    for(int i=0;i<src.N;i++){
        for(int j=0;j<i;j++){
            double t=src.element[i][j];
            src.element[i][j]=src.element[j][i];
            src.element[j][i]=t;
        }
    }
    return src;
}

//Turn a bunch of vectors into a matrix in row
matrix vectorsToMatrix(const vector *src,int num){
    matrix ret;
    ret.N=src[0].dim,ret.M=num;
    for(int i=0;i<num;i++){
        if(src[i].dim!=ret.N){
            printf("Vectors dimension mismatching!\n");
            return ret;
        }
        for(int j=0;j<src[i].dim;j++){
            ret.element[i][j]=src[i].cod[j];
        }
    }
    return ret;
}

//Turn a matrix into a bunch of vectors by rows
vector* matrixToVectors(const matrix src){
    vector* ret=(vector*)malloc(sizeof(vector)*src.M);
    for(int i=0;i<src.M;i++){
        ret[i].dim=src.N;
        for(int j=0;j<src.N;j++){
            ret[i].cod[j]=src.element[j][i];
        }
    }
    return ret;
}

//Horizontally concating two matrices src and dst into a new matrix (src,dst)
matrix matrixHorizontallyConcating(matrix src,matrix dst){
    if(src.N!=dst.N){
        printf("Dimension mismatching!\n");
        return src;
    }

    matrix ret;
    ret.N=src.N,ret.M=src.M+dst.M;
    for(int i=0;i<src.M;i++){
        for(int j=0;j<src.N;j++){
            ret.element[j][i]=src.element[j][i];
        }
    }
    for(int i=src.M;i<ret.M;i++){
        for(int j=0;j<dst.N;j++){
            ret.element[j][i]=dst.element[j][i-src.M];
        }
    }
    return ret;
}

//Vertically concating two matrices src and dst into a new matrix (src)
//                                                                (dst)
matrix matrixVerticallyConcating(matrix src,matrix dst){
    if(src.M!=dst.M){
        printf("Dimension mismatching!\n");
        return src;
    }

    matrix ret;
    ret.M=src.M,ret.N=src.N+dst.N;
    for(int i=0;i<src.M;i++){
        for(int j=0;j<src.N;j++){
            ret.element[j][i]=src.element[j][i];
        }
    }
    for(int i=0;i<ret.M;i++){
        for(int j=src.N;j<ret.N;j++){
            ret.element[j][i]=dst.element[j-src.N][i];
        }
    }
    return ret;
}

//Perform the gaussian elimination on the matrix src
//Transform the matix src into an echelon matrix
//aug: value can be only 0 or 1(1 means src is an augmented matrix;vise versa)
matrix gaussianElimination(matrix src,int aug){
    int mainVar[MAX_DIM];//maintain the corresponding line of each main variable
    memset(mainVar,-1,sizeof(mainVar));
    //cnt: the id of current line
    //i: the subscription of current unknown variable
    for(int i=0,cnt=0;i<src.M-aug&&cnt<src.N;i++){
        if(iszero(src.element[cnt][i])){
            int flag=0;
            for(int j=cnt+1;j<src.N;j++){
                if(!iszero(src.element[j][i])){
                    for(int k=0;k<src.M;k++){
                        double t=src.element[cnt][k];
                        src.element[cnt][k]=src.element[j][k];
                        src.element[j][k]=t;
                    }
                    flag=1;
                    break;
                }
            }
            if(!flag)continue;
        }

        mainVar[i]=cnt;
        double div=src.element[cnt][i];
        for(int j=i;j<src.M;j++){
            src.element[cnt][j]/=div;
        }

        for(int j=cnt+1;j<src.N;j++){
            if(iszero(src.element[j][i]))continue;
            double div=src.element[j][i];
            for(int k=i;k<src.M;k++){
                src.element[j][k]-=src.element[cnt][k]*div;
            }
        }
        cnt++;
    }
    //now the matrix is an echelon matrix

    //turn into a simplified echelon matrix
    for(int i=src.M-aug-1;i>=0;i--){
        if(mainVar[i]!=-1){
            for(int j=mainVar[i]-1;j>=0;j--){
                if(!iszero(src.element[j][i])){
                    double div=src.element[j][i];
                    for(int k=i;k<src.M;k++){
                        src.element[j][k]-=div*src.element[mainVar[i]][k];
                    }
                }
            }
        }
    }

    return src;
}

//Calculate the rank of a matrix or a bunch of vectors
int getRank(matrix src){
    src=gaussianElimination(src,0);
    int ret;
    for(int i=0;i<src.N;i++){
        int flag=1;
        for(int j=0;j<src.M;j++){
            if(src.element[i][j]){
                flag=0;
                break;
            }
        }
        if(flag)ret++;
    }
    return ret;
}

//Get inverse of matrix src
matrix getInverseMatrix(matrix src){

}

//output the vector src
void displayVector(vector src){
    printf("Dimension: %d\n",src.dim);
    for(int i=0;i<src.dim;i++){
        printf("%.2lf ",src.cod[i]);
    }
    printf("\n");
}

//output the matrix src
void displayMatrix(matrix src){
    printf("Dimension: %d x %d\n",src.N,src.M);
    for(int i=0;i<src.N;i++){
        for(int j=0;j<src.M;j++){
            printf("%.2lf ",src.element[i][j]);
        }
        printf("\n");
    }
}

//input a vector variable
void readVector(vector *src){
    printf("Input dimension of the vector:\n");
    scanf("%d",&(src->dim));
    printf("Input coordinates of the vector under natural basis:\n");
    for(int i=0;i<src->dim;i++){
        scanf("%lf",&(src->cod[i]));
    }
}

//input a vector variable
void readMatrix(matrix *src){
    printf("Input dimension of the matrix(format:N M):\n");
    scanf("%d%d",&(src->N),&(src->M));
    printf("Input each element of the matrix:\n");
    for(int i=0;i<src->N;i++){
        for(int j=0;j<src->M;j++){
            scanf("%lf",&(src->element[i][j]));
        }
    }
}