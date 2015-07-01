// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef GET_INIT_FIXED_FREE_INDEP_H
#define GET_INIT_FIXED_FREE_INDEP_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Determine the set of initially fixed and free independent variables.
 *
 * This function will drop all but one of the fixed independent variables
 * whos Lagrange multiplier is above a predefined value.
 *
 * ToDo: Finish documentation.
 *
 * @param  n         [in] Total number of variables.
 * @param  r         [in] Number of decomposed constraints.
 * @param  nu_indep  [in] Sparse vector (size == n-r) of Lagrange multipliers for
 *                   the independent variables.
 * @param  super_basic_mult_drop_tol
 *                   [in] Tolerance of nu_indep(i)/||nu_indep||inf below which
 *                   active variables will not be dropped from the superbasis.
 * @param  olevel    [in] Printing level.
 * @param  out       [out] Stream the output is printed to based on olevel.
 * @param  n_pz_X    [out] Number of dropped super basic variables (n_pz_X <= nu_indep.nz()).
 * @param  n_pz_R    [out] Number of free super basic variables (n-r == n_pz_R + n_pz_X)
 * @param  i_x_free  [out] Array (size >= n_pz_R) of indices the free independent (superbasic) variables.
 * @param  i_x_fixed [out] Array (size >= n_pz_X) of indices the droped (nonbasic) variables.
 * @param  bnd_fixed [out] Array (size >= n_pz_X) of the bounds for the dropped (nonbasic) variables.
 */
void get_init_fixed_free_indep(
  const size_type                        n
  ,const size_type                       r
  ,const SpVectorSlice                   &nu_indep
  ,const value_type                      super_basic_mult_drop_tol
  ,EJournalOutputLevel                   olevel
  ,std::ostream                          &out
  ,size_type                             *n_pz_X
  ,size_type                             *n_pz_R
  ,size_type                             i_x_free[]
  ,size_type                             i_x_fixed[]
  ,ConstrainedOptPack::EBounds  bnd_fixed[]
  );

} // end namespace MoochoPack

#endif // GET_INIT_FIXED_FREE_INDEP_H
