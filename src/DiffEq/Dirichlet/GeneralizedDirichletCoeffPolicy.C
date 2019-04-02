// *****************************************************************************
/*!
  \file      src/DiffEq/Dirichlet/GeneralizedDirichletCoeffPolicy.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Generalized Dirichlet coefficients policies
  \details   This file defines coefficients policy classes for the Generalized
             Dirichlet SDE, defined in DiffEq/Dirichlet/GeneralizedDirichlet.h.
             For general requirements on Generalized Dirichlet SDE coefficients
             policy classes see the header file.
*/
// *****************************************************************************

#include "GeneralizedDirichletCoeffPolicy.h"

walker::GeneralizedDirichletCoeffConst::GeneralizedDirichletCoeffConst(
  tk::ctr::ncomp_type ncomp,
  const std::vector< kw::sde_b::info::expect::type >& b_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappa::info::expect::type >& k_,
  const std::vector< kw::sde_c::info::expect::type >& c_,
  std::vector< kw::sde_b::info::expect::type >& b,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_c::info::expect::type >& c )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] b_ Vector used to initialize coefficient vector b
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] k_ Vector used to initialize coefficient vector k
//! \param[in] c_ Vector used to initialize coefficient vector c
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
//! \param[in,out] c Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( b_.size() == ncomp,
          "Wrong number of generalized Dirichlet SDE parameters 'b'");
  ErrChk( S_.size() == ncomp,
          "Wrong number of generalized Dirichlet SDE parameters 'S'");
  ErrChk( k_.size() == ncomp,
          "Wrong number of generalized Dirichlet SDE parameters 'kappa'");
  ErrChk( c_.size() == ncomp*(ncomp-1)/2,
          "Wrong number of generalized Dirichlet SDE parameters 'c'");

  b = b_;
  S = S_;
  k = k_;
  c = c_;
}
