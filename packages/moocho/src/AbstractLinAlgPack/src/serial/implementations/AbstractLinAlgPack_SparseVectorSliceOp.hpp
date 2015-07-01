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

#ifndef SPARSE_VECTOR_SLICE_OP_H
#define SPARSE_VECTOR_SLICE_OP_H

#include "AbstractLinAlgPack_SparseVectorOp.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"

namespace AbstractLinAlgPack {

/** @name Inline functions to call linear algebra operations for SparseVectoSlice<> objects.
  *
  * These inline functions use the same exact syntax as the DenseLinAlgPack linear algebra
  * operations.  They call the generic templated linear algebra operations for the
  * SparseVectorTemplateInterface interface.  With these functions the abreaviation
  * for the sparse vector #SV# is replaced with the standard #V#.  The
  * assignmet operation #V_SV# is replaced with the more uniform #assign#
  * identifier so as to be compatable with DenseLinAlgPack naming.  Also, the suffixes
  * for #dot#, #norm_1#, and the other level one BLAS operations is dropped.
  *
  * In order to use SparseVector<> objects with these functions the client must
  * perform an explicit conversion.  For example the following code will not
  * compile (can not perform conversion from const SparseVector<> to const SparseVectorSlice<>):
  *
  *	SparseVector<element_type> sv;
  *	// ...
  * cout << norm_1(sv);
  *
  * Instead you must perform an explicit conversion for example by using operator()() as follows:
  *
  * cout << norm_1(sv());
  */
//@{ begin Sparse BLAS operations


/// result = dot(vs_rhs1,sv_rhs2) (BLAS xDOT)
template<class T_Ele>
inline
value_type dot(const DVectorSlice& vs_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2)
{
  return dot_V_SV(vs_rhs1, sv_rhs2);
}

/// result = dot(sv_rhs1,vs_rhs2) (BLAS xDOT)
template<class T_Ele>
inline
value_type dot(const SparseVectorSlice<T_Ele>& sv_rhs1, const DVectorSlice& vs_rhs2)
{
  return dot_SV_V(sv_rhs1, vs_rhs2);
}

/// result = ||sv_rhs||1 (BLAS xASUM)
template<class T_Ele>
inline
value_type norm_1(const SparseVectorSlice<T_Ele>& sv_rhs)
{
  return norm_1_SV(sv_rhs);
}

/// result = ||sv_rhs||2 (BLAS xNRM2)
template<class T_Ele>
inline
value_type norm_2(const SparseVectorSlice<T_Ele>& sv_rhs)
{
  return norm_2_SV(sv_rhs);
}

/// result = ||sv_rhs||inf (BLAS IxAMAX)
template<class T_Ele>
inline
value_type norm_inf(const SparseVectorSlice<T_Ele>& sv_rhs)
{
  return norm_inf_SV(sv_rhs);
}

/// result = max(sv_rhs)
template<class T_Ele>
inline
value_type max(const SparseVectorSlice<T_Ele>& sv_rhs)
{
  return max_SV(sv_rhs);
}

/// result = min(sv_rhs)
template<class T_Ele>
inline
value_type min(const SparseVectorSlice<T_Ele>& sv_rhs)
{
  return min_SV(sv_rhs);
}

/// vs_lhs += alpha * sv_rhs (BLAS xAXPY)
template<class T_Ele>
inline
void Vp_StV(DVectorSlice* vs_lhs, value_type alpha, const SparseVectorSlice<T_Ele>& sv_rhs)
{
  Vp_StSV(vs_lhs, alpha, sv_rhs);
}

/// vs_lhs += alpha * op(gms_rhs1) * sv_rhs2 (BLAS xGEMV)
template<class T_Ele>
inline
void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2)
{
  Vp_StMtSV(vs_lhs, alpha, gms_rhs1, trans_rhs1, sv_rhs2);
}

/// vs_lhs += alpha * op(tri_gms_rhs1) * sv_rhs2 (BLAS xTRMV)
template<class T_Ele>
inline
void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSliceTri& tri_gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2)
{
  Vp_StMtSV(vs_lhs, alpha, tri_gms_rhs1, trans_rhs1, sv_rhs2);
}

/// vs_lhs += alpha * op(sym_gms_rhs1) * sv_rhs2 (BLAS xSYMV)
template<class T_Ele>
inline
void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSliceSym& sym_gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2)
{
  Vp_StMtSV(vs_lhs, alpha, sym_gms_rhs1, trans_rhs1, sv_rhs2);
}

/** \brief vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs
  *
  * Calls: #Vp_StMtSV(vs_lhs,alpha,M_rhs1,trans_rhs1,V_rhs2);#
  *
  * Needed for LinAlgOpPack template functions and provides
  * implementations for gms, tri, and sym matrices.
  *
  * The name Vp_StMtV could not be used because of a conflict
  * with another functions.
  */
template<class M, class T_Ele>
inline
void Vp_StMtSVS(DVectorSlice* vs_lhs, value_type alpha, const M& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2
  , value_type beta)
{
  using DenseLinAlgPack::Vp_MtV_assert_sizes;
  using DenseLinAlgPack::Vt_S;
  Vp_MtV_assert_sizes(vs_lhs->dim(), M_rhs1.rows(), M_rhs1.cols(), trans_rhs1, sv_rhs2.dim());
  if(beta == 0.0)
    *vs_lhs = 0.0;
  else if(beta != 1.0)
    Vt_S(vs_lhs,beta);
  Vp_StMtV(vs_lhs,alpha,M_rhs1,trans_rhs1,sv_rhs2);
}

/// vs_lhs = alpha * op(gms_rhs1) * sv_rhs2 + beta * vs_lhs
template<class T_Ele>
inline
void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2
  , value_type beta)
{	Vp_StMtSVS(vs_lhs, alpha, gms_rhs1, trans_rhs1, sv_rhs2, beta); }

/// vs_lhs = alpha * op(tri_gms_rhs1) * sv_rhs2 + beta * vs_lhs
template<class T_Ele>
inline
void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSliceTri& tri_gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2
  , value_type beta)
{	Vp_StMtSVS(vs_lhs, alpha, tri_gms_rhs1, trans_rhs1, sv_rhs2, beta); }

/// vs_lhs = alpha * op(sym_gms_rhs1) * sv_rhs2 + beta * vs_lhs
template<class T_Ele>
inline
void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSliceSym& sym_gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SparseVectorSlice<T_Ele>& sv_rhs2
  , value_type beta)
{	Vp_StMtSVS(vs_lhs, alpha, sym_gms_rhs1, trans_rhs1, sv_rhs2, beta); }

//@} end Sparse BLAS Operations

} // end namespace AbstractLinAlgPack

#endif // SPARSE_VECTOR_SLICE_OP_H
