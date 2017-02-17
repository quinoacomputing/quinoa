/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/*
 * This code is a wrapper layer around AZTEC.  It takes the real and imaginary
 * parts of a complex valued linear system and forms an equivalent real system,
 * calls AZTEC to solve the real system and returns the solution.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <std::string.h>
#include "az_aztec.h"
#include "az_ifpack.h"

void AZ_ifpack_iterate(double *x, double *b,
			int *options, double *params,
			double *status, int *proc_config,
			AZ_MATRIX *Amat )

/*******************************************************************************

  Author:          Mike Heroux, SNL, 9222
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the initial guess, On output 
                   contains the solution to 
                   the linear system.

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

 Amat:             The matrix operator, stored as an AZ_MATRIX structure.

  Internal Parameter list:
  ========================

 x:                Komplex version of initial guess and solution.
 b:                Komplex version of RHS.
 Prec:             Preconditioner stored as an AZ_PRECOND structure.

  Overview
  ========

*******************************************************************************/


{

AZ_PRECOND *Prec;    /* Structure representing entire preconditioner.        */
                    /*                                                      */

   AZ_ifpack_prec_create (x, b, options,  params, proc_config, Amat, &Prec);

   /* solve linear system using Aztec. */

   AZ_iterate(x, b, options, params, status, proc_config, Amat, Prec, NULL);

   AZ_ifpack_prec_destroy (options,  params, proc_config, Amat, Prec);

/* AZ_ifpack_iterate*/
}
