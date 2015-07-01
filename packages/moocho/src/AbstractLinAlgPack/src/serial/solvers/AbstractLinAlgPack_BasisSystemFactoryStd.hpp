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

#ifndef SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H
#define SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_BasisSystemFactory.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace AbstractLinAlgPack {

/** \brief Default implementation for <tt>BasisSystemPermDirectSparse</tt> obejcts
 * using <tt>DirectSparseSolver</tt> object.
 *
 * Several direct sparse solvers are supported by default.  These include:
 * <ul>
 * <li> DENSE (using LAPACK xGETRF())
 * <li> MA28
 * <li> MA48 (using MA28 for BasisSystemPerm::select_basis()) (not yet)
 * <li> SuperLU
 * </ul>
 *
 * These solvers are supported only if the proper macros are defined.
 * 
 * ToDo: Create a DirectSparseSolverFactory interface and use this
 * to allow clients to add new DirectSparseSolvers ...
 *
 */
class BasisSystemFactoryStd
  : public AbstractLinAlgPack::BasisSystemFactory
{
public:

  /** \brief . */
  BasisSystemFactoryStd(); // ToDo: Add arguments!

  /** @name Overridden from BasisSystemFactory */
  //@{

  /** \brief . */
  void set_options( const options_ptr_t& options );
  /** \brief . */
  const options_ptr_t& get_options() const;

  //@}

  /** @name Overridden from AbstractFactory */
  //@{

  /** \brief . */
  obj_ptr_t create() const;

  //@}

private:

  // ////////////////////////
  // Private types

  enum EDirectLinearSolverType { LA_DENSE, LA_MA28, LA_MA48, LA_SUPERLU };

  // ////////////////////////
  // Private data members

  mutable EDirectLinearSolverType  direct_linear_solver_type_;
  options_ptr_t                    options_;

  // ////////////////////////
  // Private member functions
  
  void read_options() const;

}; // end class BasisSystemFactoryStd

}  // end namespace AbstractLinAlgPack

#endif // SPARSE_SOLVER_PACK_BASIS_SYSTEM_FACTORY_H
