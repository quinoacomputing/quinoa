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

#include "ConstrainedOptPack_MatrixIdentConcat.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_AssertOp.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp" // RAB: 12/20/2002: Must be after DenseLinAlgPack_LinAlgOpPack.hpp
                                                                        // for Intel C++ 5.0 (Windows 2000).  Namespaces
                                                                        // and lookup rules don't work properly with this compiler.
#include "Teuchos_FancyOStream.hpp"

namespace {

// Get a view of a vector (two versions)

inline
Teuchos::RCP<const AbstractLinAlgPack::Vector>
get_view(
  const AbstractLinAlgPack::Vector  &v
  ,const RangePack::Range1D         &rng
  )
{
   return v.sub_view(rng);
}

inline
Teuchos::RCP<const AbstractLinAlgPack::SpVectorSlice>
get_view(
  const AbstractLinAlgPack::SpVectorSlice &v
  ,const RangePack::Range1D               &rng
  )
{
  return Teuchos::rcp( new AbstractLinAlgPack::SpVectorSlice( v(rng) ) );
}

// Matrix-vector multiplication

template<class V>
void mat_vec(
  AbstractLinAlgPack::VectorMutable        *y
  ,AbstractLinAlgPack::value_type          a
  ,AbstractLinAlgPack::value_type          alpha
  ,const AbstractLinAlgPack::MatrixOp      &D
  ,BLAS_Cpp::Transp                        D_trans
  ,const DenseLinAlgPack::Range1D          &D_rng
  ,const DenseLinAlgPack::Range1D          &I_rng
  ,BLAS_Cpp::Transp                        M_trans
  ,const V                                 &x
  ,AbstractLinAlgPack::value_type          b
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;

  //
  // M = [ alpha*op(D) ] or [      I       ]
  //     [     I       ]    [ alpha*op(D)  ]
  //
  if( M_trans == no_trans ) {
    // 
    // y = b*y + a*M*x
    // =>
    // y(D_rng) = b*y(D_rng) + a*alpha*op(D)*x
    // y(I_rng) = b*y(I_rng) + a*x
    //
    AbstractLinAlgPack::VectorSpace::vec_mut_ptr_t
      y_D = y->sub_view(D_rng),
      y_I = y->sub_view(I_rng);
    // y(D_rng) = b*y(D_rng) + a*alpha*op(D)*x
    AbstractLinAlgPack::Vp_StMtV( y_D.get(), a*alpha, D, D_trans, x, b );
    // y(I_rng) = b*y(I_rng) + a*x
    LinAlgOpPack::Vt_S( y_I.get(), b );
    LinAlgOpPack::Vp_StV( y_I.get(), a, x );
  }
  else {
    //
    // y = b*y + a*M'*x
    // =>
    // y = b*y + a*alpha*op(D')*x(D_rng) + a*x(I_rng)
    //
    LinAlgOpPack::Vp_StMtV( y, a*alpha, D, trans_not(D_trans), *get_view(x,D_rng), b );
    LinAlgOpPack::Vp_StV( y, a, *get_view(x,I_rng) );
  }
}

} // end namespace

namespace ConstrainedOptPack {

// Overridden from MatrixBase

size_type MatrixIdentConcat::rows() const
{
  return this->D_rng().size() + this->I_rng().size();
}

size_type MatrixIdentConcat::cols() const
{
  return this->I_rng().size();
}

size_type MatrixIdentConcat::nz() const
{
  const MatrixOp& D = this->D();
  return D.nz() + BLAS_Cpp::cols(D.rows(),D.cols(),D_trans()); // D and I
}

// Overridden from MatrixOp

std::ostream& MatrixIdentConcat::output(std::ostream& out_arg) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  const Range1D           D_rng   = this->D_rng();
  const BLAS_Cpp::Transp  D_trans = this->D_trans();
  *out << "Converted to dense =\n";
  MatrixOp::output(*out);
  *out << "This is a " << this->rows() << " x " << this->cols()
    << " general matrix / identity matrix concatenated matrix object ";
  if( D_rng.lbound() == 1 ) {
    if( D_trans == BLAS_Cpp::no_trans )
      *out << "[ alpha*D; I ]";
    else
      *out << "[ alpha*D'; I ]";
  }
  else {
    if( D_trans == BLAS_Cpp::no_trans )
      *out << "[ I; alpha*D ]";
    else
      *out << "[ I; alpha*D' ]";
  }
  *out << " where alpha and D are:";
  tab.incrTab();
  *out << "\nalpha = " << this->alpha();
  *out << "\nD =\n" << D();
  return out_arg;
}

void MatrixIdentConcat::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  , const Vector& x, value_type b
  ) const
{
  AbstractLinAlgPack::Vp_MtV_assert_compatibility(y,*this,M_trans,x);
  mat_vec( y, a, alpha(), D(), D_trans(), D_rng(), I_rng(), M_trans, x, b );
}

void MatrixIdentConcat::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b
  ) const
{

  AbstractLinAlgPack::Vp_MtV_assert_compatibility(y,*this,M_trans,x);
  mat_vec( y, a, alpha(), D(), D_trans(), D_rng(), I_rng(), M_trans, x, b );
}

} // end namespace ConstrainedOptPack
