#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>



// gcc -Wall -Werror ERAM.c -o ERAM.out -lm -lblas -llapacke



double norme(double *v, int n) {
  double res = 0;
  for(int i = 0; i < n; i++) {
    res += pow(v[i], 2);
  }
  return sqrt(res);
}

void VectScal(double *v, int n, int scal, double *res) {
  for(int i = 0; i < n; i++) {
    res[i] = v[i] * scal;
  }
}

double VectVect(double *v1, double* v2, int n) {
  double res = 0;
  for(int i = 0; i < n; i++) {
    res += v1[i] * v2[i];
  }
  return res;
}

void MatVect(double *m, double *v, int n, double *res) {
  for(int i = 0; i < n; i++) {
    res[i] = 0;

    for(int j = 0; j < n; j++) {
      res[i] += m[i + j*n] * v[j];
    }
  }
}

void VectAdd(double *v1, double *v2, int n, double *res) {
  for(int i = 0; i < n; i++) {
    res[i] = v1[i] + v2[i];
  }
}

void VectSub(double *v1, double *v2, int n, double *res) {
  for(int i = 0; i < n; i++) {
    res[i] = v1[i] - v2[i];
  }
}


double arnoldi(double *A, int n, int m, double *q, double *Q, double *H) {
  double q_norme = norme(q, n);
  
  for(int i = 0; i < n; i++) {
    Q[i] = q[i] / q_norme;
  }

  double *tmp1 = malloc(n * sizeof(double));
  double *tmp2 = malloc(n * sizeof(double));


  for(int k = 1; k <= m; k++) {
    MatVect(A, &(Q[(k-1)*n]), n, tmp1);

    for(int j = 1; j <= k; j++) {
      H[(j-1) + (k-1)*n] = VectVect(tmp1, &Q[(j-1)*n], n);
    }

    for(int j = 1; j <= k; j++) {
      VectScal(&(Q[(j-1)*n]), n, H[(j-1) + (k-1)*n], tmp2);
      VectSub(tmp1, tmp2, n, tmp1);
    }
    
    if(k == m) {
      return norme(tmp1, n);
    }

    H[k + (k-1)*n] = norme(tmp1, n);

    VectScal(tmp1, 1. / H[k + (k-1)*n], n, &(Q[k*n]));
  }

  return -1;
}

void eram(double *A, int n, int m, double *q, int k, double tol, double *wr, double* wi, double *vr) { 
  double *Q = malloc(m * n * sizeof(double));
  double *H = malloc(m * m * sizeof(double));

  //double *wr = malloc(m * sizeof(double));
  //double *wi = malloc(m * sizeof(double));
  //double *vr = malloc(m * m * sizeof(double));

  double residu = arnoldi(A, n, m, q, Q, H);
  LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', m, H, m, wr, wi, NULL, 1, vr, m);

 // trie 
 // calcul de vecteur propre


  while( residu > tol) {
    // creer le nouveau vecteur
    //copie 101->105
  }
}

int main(int argc, char **argv) {
  
  return 0;
}
