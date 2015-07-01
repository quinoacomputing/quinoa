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

#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianFull.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace ConstrainedOptPack {

void QPSchurInitKKTSystemHessianFull::initialize_kkt_system(
  const Vector    &g
  ,const MatrixOp   &G
  ,value_type           etaL
  ,const Vector   *dL
  ,const Vector   *dU
  ,const MatrixOp   *F
  ,BLAS_Cpp::Transp     trans_F
  ,const Vector   *f
  ,const Vector   *d
  ,const Vector   *nu
  ,size_type            *n_R
  ,i_x_free_t           *i_x_free
  ,i_x_fixed_t          *i_x_fixed
  ,bnd_fixed_t          *bnd_fixed
  ,j_f_decomp_t         *j_f_decomp
  ,DVector               *b_X
  ,Ko_ptr_t             *Ko
  ,DVector               *fo
  ) const
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  using LinAlgOpPack::V_mV;

  // Validate type of and convert G
  const MatrixSymOpNonsing&
    G_sym = dyn_cast<const MatrixSymOpNonsing>(G);

  const size_type nd = g.dim();
  
  // n_R
  *n_R = nd;
  // i_x_free[i-1] = i, i = 1...nd
  i_x_free->resize(0);
  // i_x_fixed[0] = nd+1
  i_x_fixed->resize(1);
  (*i_x_fixed)[0] = nd+1;
  // bnd_fixed[0] = LOWER
  bnd_fixed->resize(1);
  (*bnd_fixed)[0] = LOWER;
  // j_f_decomp[] = empty
  j_f_decomp->resize(0);
  // b_X = etaL
  b_X->resize(1);
  (*b_X)[0] = etaL;
  // Ko = G
  *Ko = Teuchos::rcp(&G_sym,false); // Not dynamically allocated so don't delete!
  // fo = -g
  V_mV(fo,VectorDenseEncap(g)());
}

} // end namesapce ConstrainedOptPack
