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

#ifndef EXAMPLE_BASIS_SYSTEM_H
#define EXAMPLE_BASIS_SYSTEM_H

#include "NLPInterfacePack_Types.hpp"
#include "AbstractLinAlgPack_BasisSystemComposite.hpp"

namespace NLPInterfacePack {

/** \brief Subclass of BasisSystem for example NLP.
 *
 * ToDo: Finish documentation!
 */
class ExampleBasisSystem
  : public AbstractLinAlgPack::BasisSystemComposite
{
public:

  /// Calls <tt>this->initialize()</tt>
  ExampleBasisSystem(
    const VectorSpace::space_ptr_t       &space_x
    ,const Range1D                       &var_dep
    ,const Range1D                       &var_indep
    );
  
  /** \brief Initialize given the vector space for the dependent and independent variables.
   *
   * @param  space_x   [in]
   * @param  var_dep   [in]
   * @param  var_indep [in]
   *
   * ToDo: Finish documentation!
   */
  void initialize(
    const VectorSpace::space_ptr_t       &space_x
    ,const Range1D                       &var_dep
    ,const Range1D                       &var_indep
    );

  /** @name Overridden from BasisSystemComposite */
  //@{

  /** \brief . */
  void update_D(
    const MatrixOpNonsing       &C
    ,const MatrixOp             &N
    ,MatrixOp                   *D
    ,EMatRelations               mat_rel
    ) const;

  //@}

private:
  // Not defined and not to be called!
  ExampleBasisSystem();

}; // end class ExampleBasisSystem

} // end namespace NLPInterfacePack

#endif // EXAMPLE_BASIS_SYSTEM_H
