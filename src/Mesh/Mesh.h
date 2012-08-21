// -----------------------------------------------------------------------------
// \file    src/Mesh/Mesh.h
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Unstructured mesh
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

void prepmsh(int *npoin, int *nbpoin, int *nelem, double **coord, int **bpg,
             int **binpoel, int **inpoel, int **esup1, int **esup2, int **psup1,
             int **psup2, int **bpsup1, int **bpsup2, int **esupel1,
             int **esupel2, int **esuel, int**bptags, int **betags, double **Ae,
             double **dNx, double **dNy, double **dete, double **sqrtAp,
             double *minsqrtAp);
