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

#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU

#ifndef SSP_SUPERLU_SOLVER_H
#define SSP_SUPERLU_SOLVER_H

#include "Teuchos_RCP.hpp"

namespace SuperLUPack {

/** \brief Abstract interface to SuperLU.
 *
 * This class is abstract so that none of the SuperLU header files need
 * be presented to the client.  SuperLU uses some very bad software engineering
 * by injecting into the global namespace lots of inappropriate names.  All
 * SuperLU global names should be prefixed by SLU_ or SuperLU_ but they are not
 */
class SuperLUSolver {
public:

  /** @name Public Types */
  //@{

  /** \brief Abstract class for objects that represent the factorization
    * structure of a particular class of matrices.
    *
    * This structure can be reused over and over again for factorizing matrices
    * with the same structure but different nonzero elements.
    */
  class FactorizationStructure {
  public:
    /** \brief . */
    virtual ~FactorizationStructure() {}
  };

  /** \brief Abstract class for objects that represent the factorization
    * nonzeos of a particular matrix.
    */
  class FactorizationNonzeros {
  public:
    /** \brief . */
    virtual ~FactorizationNonzeros() {}
  };

  //@}

  /** \brief . */
  virtual ~SuperLUSolver() {}

  /** @name Static members */
  //@{

  /** \brief . */
  static Teuchos::RCP<SuperLUSolver>                         create_solver();
  /** \brief . */
  static Teuchos::RCP<SuperLUSolver::FactorizationStructure> create_fact_struct();
  /** \brief . */
  static Teuchos::RCP<SuperLUSolver::FactorizationNonzeros>  create_fact_nonzeros();

  //@}

  /** @name Basis selection and factorization */
  //@{

  /** \brief Analyze and factor the matrix, finding a basis in the process
   *
   * param  m     [in] Number of rows in A
   * param  n     [in] Number of columns in A
   * param  nz    [in] Number of nonzeros in A
   * param  a_val [in] Array (length nz) of the values of A
   * param  a_row_i
   *              [in] Array (length nz) of the rows indexes of A (zero-based
   *              for SuperLU)
   * param  a_col_ptr
   *              [in] Array (length n+1) Column pointers for starts
   *              of columns for A in a_val[] and a_col_ptr[]
   * param  fact_struc
   *              [out] Factorization structure for basis of A
   * param  fact_nonzeros
   *              [out] Factorization nonzeros for basis of A
   * param  perm_r
   *              [out] Array (length m) of row permutations for
   *              basis selection of A (zero-based).
   * param  perm_c
   *              [out] Array (length n) of column permutations
   *              for basis selection of A (zero-based).
   * param  rank  [out] Rank of the basis of A selected.
   *
   * Preconditions:<ul>
   * <li> <tt>m >= n</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> Lots of others???
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  virtual void analyze_and_factor(
    int                         m
    ,int                        n
    ,int                        nz
    ,const double               a_val[]
    ,const int                  a_row_i[]
    ,const int                  a_col_ptr[]
    ,FactorizationStructure     *fact_struct
    ,FactorizationNonzeros      *fact_nonzeros
    ,int                        perm_r[]
    ,int                        perm_c[]
    ,int                        *rank
    ) = 0;
  
  /** \brief Refactor the same basis.
   *
   * ToDo: Finish documentation!
   */
  virtual void factor(
    int                             m
    ,int                            n
    ,int                            nz
    ,const double                   a_val[]
    ,const int                      a_row_i[]
    ,const int                      a_col_ptr[]
    ,const FactorizationStructure   &fact_struct
    ,FactorizationNonzeros          *fact_nonzeros
    ) = 0;

  //@}
  
  /** @name Solve for linear systems */
  //@{

  /** \brief Solve a set of linear systems.
   *
   * ToDo: Finish documentation!
   */
  virtual void solve(
    const FactorizationStructure    &fact_struct
    ,const FactorizationNonzeros    &fact_nonzeros
    ,bool                           transp
    ,int                            n
    ,int                            nrhs
    ,double                         rhs[]
    ,int                            ldrhs
    ) const = 0;

  //@}

}; // end class SuperLUSolver

} // end namespace SuperLUPack

#endif // SSP_SUPERLU_SOLVER_H

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
