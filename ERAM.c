#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <TError.h>



// gcc -Wall -Werror ERAM.c -o ERAM.out -lm -lblas -llapacke -std=gnu11



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

void MatVect(double *m, double *v, int w, int h, double *res) {
  for(int i = 0; i < h; i++) {
    res[i] = 0;

    for(int j = 0; j < w; j++) {
      res[i] += m[i + j*h] * v[j];
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
    MatVect(A, &(Q[(k-1)*n]), n, n, tmp1);

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

double mod(double a, double b) {
  return sqrt(pow(a,2) + pow(b, 2));
}

void customBubbleSort(double* wr, double* wi, double* vr, int m) {
  double swapwr;
  double swapwi;
  double swapvr[m];

  double max, tmp;
  int index;

  for(int i = 0; i < m - 1; i++) {
    max = mod(wr[i], wi[i]);
    index = i;
    for(int j = i+1; j < m; j++) {
      tmp = mod(wr[j], wi[j]);
      if (tmp > max) {
        max = tmp;
        index = j;
      }
    }

    if ( index != i) {
      swapwr = wr[i];
      swapwi = wi[i];
      memcpy(swapvr, &(vr[i*m]), m * sizeof(double));
      
      wr[i] = wr[index];
      wi[i] = wi[index];
      memcpy(&(vr[i*m]), &(vr[index*m]), m * sizeof(double));

      wr[index] = swapwr;
      wi[index] = swapwi;
      memcpy(&(vr[index*m]), swapvr, m * sizeof(double));
    }
  }
}

void eram(double *A, int n, int m, double *q, int k, double tol, double *wr, double* wi, double *vr) { 
  double *Q = malloc(m * n * sizeof(double));
  double *H = malloc(m * m * sizeof(double));

  double *wrt = malloc(m * sizeof(double));
  double *wit = malloc(m * sizeof(double));
  double *vrt = malloc(m * m * sizeof(double));

  if(!Q || !H || !wrt || !wit || !vrt) {
    fprintf(stderr, "AAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
    exit(1);
  }
  
  // attention dépassement mem quelque part ...
  double residu = arnoldi(A, n, m, q, Q, H);
  LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', m, H, m, wrt, wit, NULL, 1, vrt, m);
  
  //customBubbleSort(wrt, wit, vrt, m);

  // MATMAT possible (et mieux ?)
  for(int i = 0; i < m; i++) {
    MatVect(Q, &(vrt[i*m]), m, n, &(vr[i*n]));
  }


  while( residu > tol ) {
    memcpy(q, vr, n * sizeof(double));
    for(int i = 1; i < m; i++) {
      VectAdd(q, &(vr[i*n]), n, q);
    }





    residu = arnoldi(A, n, m, q, Q, H);
    LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', m, H, m, wrt, wit, NULL, 1, vrt, m);

    //customBubbleSort(wrt, wit, vrt, m);

    // MATMAT possible (et mieux ?)
    for(int i = 0; i < m; i++) {
      MatVect(Q, &(vrt[i*m]), m, n, &(vr[i*n]));
    }
  }

  free(Q);
  free(H);

  free(wrt);
  free(wit);
  free(vrt);
}

int main(int argc, char **argv) {
  TError_init();

  int n = 6, m = 3, k = m;
  double tol = 0.1;
  double *A = malloc(n * n * sizeof(double));
  double *q = malloc(n * sizeof(double));

  srand48(0);

  for (int i = 0; i < n*n; i++) {
    A[i] = floor(drand48() * 100);
    printf("%lf ", A[i]);

    if(i < n) {
      q[i] = floor(drand48() * 100);
    }
  }
  printf("\n");


  double *wr = malloc(m * sizeof(double));
  double *wi = malloc(m * sizeof(double));
  double *vr = malloc(m * n * sizeof(double));
  
  eram(A, n, m, q, k, tol, wr, wi, vr);

  for(int i = 0; i < m; i++) {
    printf("\n%lf + %lfi\n\t", wr[i], wi[i]);
    for(int j = 0; j < n; j++) {
      printf("%lf ", vr[i*n + j]);
    }
  }

  return 0;
}
