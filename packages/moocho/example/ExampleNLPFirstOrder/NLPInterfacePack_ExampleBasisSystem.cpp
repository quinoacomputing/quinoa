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

#include "NLPInterfacePack_ExampleBasisSystem.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_MatrixComposite.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

namespace NLPInterfacePack {
 
ExampleBasisSystem::ExampleBasisSystem(
  const VectorSpace::space_ptr_t       &space_x
  ,const Range1D                       &var_dep
  ,const Range1D                       &var_indep
  )
  :BasisSystemComposite(
    space_x
    ,var_dep
    ,var_indep
    ,space_x->sub_space(var_dep)
    ,Teuchos::rcp(
      new Teuchos::AbstractFactoryStd<MatrixOpNonsing,MatrixSymDiagStd>())       // C
    ,Teuchos::rcp(
      new Teuchos::AbstractFactoryStd<MatrixSymOp,MatrixSymDiagStd>())           // D'*D
    ,Teuchos::rcp(
      new Teuchos::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymDiagStd>())    // S
    ,Teuchos::rcp(
      new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSymDiagStd>())              // D
    )
{}
  
void ExampleBasisSystem::initialize(
  const VectorSpace::space_ptr_t       &space_x
  ,const Range1D                       &var_dep
  ,const Range1D                       &var_indep
  )
{
  namespace mmp = MemMngPack;
  TEUCHOS_TEST_FOR_EXCEPTION(
    space_x.get() == NULL, std::invalid_argument
    ,"ExampleBasisSystem::initialize(...) : Error, space_x must be specified!"
    );
  BasisSystemComposite::initialize(
    space_x
    ,var_dep
    ,var_indep
    ,space_x->sub_space(var_dep)
    ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOpNonsing,MatrixSymDiagStd>())      // C
    ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOp,MatrixSymDiagStd>())          // D'*D
    ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymDiagStd>())   // S
    ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSymDiagStd>())             // D
    );
}

void ExampleBasisSystem::update_D(
  const MatrixOpNonsing   &C
  ,const MatrixOp         &N
  ,MatrixOp               *D
  ,EMatRelations          mat_rel
  ) const
{
  using Teuchos::dyn_cast;

  TEUCHOS_TEST_FOR_EXCEPTION(
    D == NULL, std::logic_error
    ,"ExampleBasisSystem::update_D(...): Error!" );

  const MatrixSymDiagStd
    &C_aggr = dyn_cast<const MatrixSymDiagStd>(C),
    &N_aggr = dyn_cast<const MatrixSymDiagStd>(N);
  MatrixSymDiagStd
    &D_sym_diag = dyn_cast<MatrixSymDiagStd>(*D);
  if( D_sym_diag.rows() != C.rows() )
    D_sym_diag.initialize(
      this->space_x()->sub_space(this->var_dep())->create_member()
      );
  AbstractLinAlgPack::ele_wise_divide(                           // D_diag = - N_diag ./ C_diag
    -1.0, N_aggr.diag(), C_aggr.diag(), &D_sym_diag.diag() );  // ...
}

} // end namespace NLPInterfacePack
