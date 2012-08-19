// -----------------------------------------------------------------------------
// \file    src/LinearAlgebra/CG.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Conjugate gradient solver
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

int cg_solve( sparsemat *A, double *b, double *x, int use_pc );
