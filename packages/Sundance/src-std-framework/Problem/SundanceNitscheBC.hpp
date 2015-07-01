/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_NITSCHEBC_H
#define SUNDANCE_NITSCHEBC_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"


namespace Sundance
{
class CellFilter;
class QuadratureFamily;

/**
 * This function forms the expressions that apply the Dirichlet 
 * BC \f$u=u_{BC}\f$ via Nitsche's method for the Poisson operator.
 */
Expr NitschePoissonDirichletBC(int dim,
  const CellFilter& cells,
  const QuadratureFamily& quad,
  const Expr& kappa,
  const Expr& v,
  const Expr& u,
  const Expr& uBC,
  const double& gamma);

/**
 * This function forms the expressions that apply the Dirichlet 
 * BC \f${\bf u}={\bf u}_{BC}\f$ via Nitsche's method for the Stokes operator.
 * 
 * \param cells the surface on which the BC is to be applied
 * \param quad the quadrature rule to be used
 * \param nu the viscosity, which may be a function of velocity
 * \param v velocity test function
 * \param q pressure test function
 * \param u velocity unknown function
 * \param p pressure unknown function
 * \param uBC specified velocity at surface
 */
Expr NitscheStokesNoSlipBC(const CellFilter& cells,
  const QuadratureFamily& quad,
  const Expr& nu,
  const Expr& v,
  const Expr& q,
  const Expr& u,
  const Expr& p,
  const Expr& uBC,
  const double& gamma1,
  const double& gamma2
  );




}


#endif
