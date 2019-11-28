#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

int main(int argc, char ** argv) {
  int m = 3;
  double H[] = 
  {1., 0., 0.,
   0., 2., 0.,
   0., 0., 3.};
  
  double wr[m], wi[m], vr[m*m];

  LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', m, H, m, wr, wi, NULL, 1, vr, m);

  for(int i = 0; i < m; i++) {
    printf("%lf + %lfi\n[ ", wr[i], wi[i]);
    for(int j = 0; j < m; j++) {
      printf("%lf , ", vr[i*m+j]);
    }
    printf("]\n-----------------\n");
  }
}
