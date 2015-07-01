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
//

#include "AbstractLinAlgPack_MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

MatrixOpNonsingThyra::MatrixOpNonsingThyra()
{}

MatrixOpNonsingThyra::MatrixOpNonsingThyra(
  const Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
  ,BLAS_Cpp::Transp                                                            thyra_linear_op_trans
  )
{
  this->initialize(thyra_linear_op_ns,thyra_linear_op_trans);
}

void MatrixOpNonsingThyra::initialize(
  const Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
  ,BLAS_Cpp::Transp                                                            thyra_linear_op_trans
  )
{
  namespace mmp = MemMngPack;
  TEUCHOS_TEST_FOR_EXCEPTION(
    thyra_linear_op_ns.get()==NULL, std::invalid_argument
    ,"MatrixOpNonsingThyra::initialize(thyra_linear_op_ns): Error!"
    );
  MatrixOpThyra::initialize(thyra_linear_op_ns,thyra_linear_op_trans);
}

Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> > 
MatrixOpNonsingThyra::set_uninitialized()
{
  Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >
    tmp_thyra_linear_op_ns = thyra_linear_op_ns();
  MatrixOpThyra::set_uninitialized();
  return tmp_thyra_linear_op_ns;
}

Teuchos::RCP<const Thyra::LinearOpWithSolveBase<value_type> >
MatrixOpNonsingThyra::thyra_linear_op_ns() const
{
  return Teuchos::rcp_dynamic_cast<const Thyra::LinearOpWithSolveBase<value_type> >(this->thyra_linear_op());
}

// Overridden from MatrixOp (needed to remove ambiguities)

MatrixOp::mat_mut_ptr_t
MatrixOpNonsingThyra::clone()
{
  return this->MatrixOpThyra::clone();
}

// Overridden from MatrixNonsing

void MatrixOpNonsingThyra::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  ,const Vector& v_rhs2
  ) const
{
  using Teuchos::dyn_cast;
  using BLAS_Cpp::trans_trans;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    v_lhs==NULL, std::invalid_argument
    ,"MatrixOpThyra::Vp_StMtV(...): Error!"
    );
#endif
  *v_lhs = 0.0; // Must initialize before sending to solve(...)!
  VectorMutableThyra &v_thyra_lhs = dyn_cast<VectorMutableThyra>(*v_lhs);
  Teuchos::RCP<Thyra::VectorBase<value_type> > thyra_vec_lhs = v_thyra_lhs.set_uninitialized();
  Thyra::solve<value_type>(
    *thyra_linear_op_ns()
    ,trans_trans(trans_rhs1,thyra_linear_op_trans())==BLAS_Cpp::no_trans ? Thyra::NOTRANS : Thyra::TRANS  // M_trans
    ,*dyn_cast<const VectorMutableThyra>(v_rhs2).thyra_vec()                                              // y
    ,thyra_vec_lhs.ptr()                                                                                  // x
    );
  v_thyra_lhs.initialize(thyra_vec_lhs);
}

void MatrixOpNonsingThyra::M_StInvMtM(
  MatrixOp* m_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  ) const
{
  MatrixNonsing::M_StInvMtM(m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2); // ToDo: Specialize!
}

// Overridden from MatrixOpNonsing

MatrixOpNonsing::mat_mwons_ptr_t
MatrixOpNonsingThyra::clone_mwons() const
{
  return Teuchos::null; // ToDo: Add a clone function to Thyra::LinearOpWithSolveBase???
  //return Teuchos::rcp(new MatrixOpNonsingThyra(thyra_linear_op_ns()->clone_lows()));
}

} // end namespace AbstractLinAlgPack
