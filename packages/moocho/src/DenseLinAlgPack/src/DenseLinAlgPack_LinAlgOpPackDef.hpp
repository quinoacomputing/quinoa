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

#ifndef LIN_ALG_OP_PACK_DEF_H
#define LIN_ALG_OP_PACK_DEF_H

#include "DenseLinAlgPack_LinAlgOpPackDecl.hpp"	// also includes some inline function definitions
#include "DenseLinAlgPack_AssertOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace LinAlgOpPack {

using BLAS_Cpp::rows;
using BLAS_Cpp::cols;

// Inject assert functions
using DenseLinAlgPack::assert_gms_lhs;
using DenseLinAlgPack::Vp_V_assert_sizes;
using DenseLinAlgPack::VopV_assert_sizes;
using DenseLinAlgPack::Mp_M_assert_sizes;
using DenseLinAlgPack::MopM_assert_sizes;
using DenseLinAlgPack::Vp_MtV_assert_sizes;
using DenseLinAlgPack::MtV_assert_sizes;
using DenseLinAlgPack::MtM_assert_sizes;

// Inject names of base linear algebra functions for DenseLinAlgPack.
// Note that this is neccesary in MS VC++ 5.0 because
// it does not perform name lookups properly but it
// is not adverse to the standard so it is a portable
// fix.
using DenseLinAlgPack::assign;
using DenseLinAlgPack::Vt_S;
using DenseLinAlgPack::Vp_StV;
using DenseLinAlgPack::Vp_StMtV;
using DenseLinAlgPack::Mt_S;
using DenseLinAlgPack::Mp_StM;
using DenseLinAlgPack::Mp_StMtM;

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Vectors

// //////////////////////////////////////////////////////////////////////////////
// += operations 

// //////////////////////////////////////////////////////////////////////////////
// operations with DVector as lhs

// v_lhs = V_rhs.
template <class V>
void assign(DVector* v_lhs, const V& V_rhs) {
  v_lhs->resize(V_rhs.dim());
  (*v_lhs) = 0.0;
  Vp_V(&(*v_lhs)(),V_rhs);
}

// v_lhs = alpha * V_rhs.
template <class V>
void V_StV(DVector* v_lhs, value_type alpha, const V& V_rhs) {
  v_lhs->resize(V_rhs.dim());
  (*v_lhs) = 0.0;
  Vp_StV(&(*v_lhs)(),alpha,V_rhs);
}

// v_lhs = V1_rhs1 + V2_rhs2.
template <class V1, class V2>
void V_VpV(DVector* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
  VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
  v_lhs->resize(V1_rhs1.dim());
  (*v_lhs) = 0.0;
  DVectorSlice vs_lhs(*v_lhs);
  Vp_V(&vs_lhs,V1_rhs1);
  Vp_V(&vs_lhs,V2_rhs2);
}


// v_lhs = V_rhs1 - V_rhs2.
template <class V1, class V2>
void V_VmV(DVector* v_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
  VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
  v_lhs->resize(V1_rhs1.dim());
  (*v_lhs) = 0.0;
  DVectorSlice vs_lhs(*v_lhs);
  Vp_V(&vs_lhs,V1_rhs1);
  Vp_StV(&vs_lhs,-1.0,V2_rhs2);
}


// v_lhs = alpha * V_rhs1 + vs_rhs2.
template <class V>
void V_StVpV(DVector* v_lhs, value_type alpha, const V& V_rhs1
  , const DVectorSlice& vs_rhs2)
{
  VopV_assert_sizes(V_rhs1.dim(),vs_rhs2.dim());
  (*v_lhs) = vs_rhs2;
  Vp_StV(&(*v_lhs)(),alpha,V_rhs1);
}

// ///////////////////////////////////////////////////////////////////////////
// operations with DVectorSlice as lhs

// vs_lhs = V_rhs.
template <class V>
void assign(DVectorSlice* vs_lhs, const V& V_rhs) {
  Vp_V_assert_sizes( vs_lhs->dim(), V_rhs.dim() );
  (*vs_lhs) = 0.0;
  Vp_V(vs_lhs,V_rhs);
}

// vs_lhs = alpha * V_rhs.
template <class V>
void V_StV(DVectorSlice* vs_lhs, value_type alpha, const V& V_rhs) {
  Vp_V_assert_sizes( vs_lhs->dim(), V_rhs.dim() );
  (*vs_lhs) = 0.0;
  Vp_StV(vs_lhs,alpha,V_rhs);
}

// vs_lhs = V1_rhs1 + V2_rhs2.
template <class V1, class V2>
void V_VpV(DVectorSlice* vs_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
  VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
  Vp_V_assert_sizes( vs_lhs->dim(), V1_rhs1.dim() );
  (*vs_lhs) = 0.0;
  Vp_V(vs_lhs,V1_rhs1);
  Vp_V(vs_lhs,V2_rhs2);
}

// vs_lhs = V_rhs1 - V_rhs2.
template <class V1, class V2>
void V_VmV(DVectorSlice* vs_lhs, const V1& V1_rhs1, const V2& V2_rhs2) {
  VopV_assert_sizes(V1_rhs1.dim(),V2_rhs2.dim());
  Vp_V_assert_sizes( vs_lhs->dim(), V1_rhs1.dim() );
  (*vs_lhs) = 0.0;
  Vp_V(vs_lhs,V1_rhs1);
  Vp_StV(vs_lhs,-1.0,V2_rhs2);
}

// vs_lhs = alpha * V_rhs1 + vs_rhs2.
template <class V>
void V_StVpV(DVectorSlice* vs_lhs, value_type alpha, const V& V_rhs1
  , const DVectorSlice& vs_rhs2)
{
  VopV_assert_sizes(V_rhs1.dim(),vs_rhs2.dim());
  (*vs_lhs) = vs_rhs2;
  Vp_StV(vs_lhs,alpha,V_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// Level 1 BLAS for Matrices

// //////////////////////////////////////////////////////////////////////////////
// += operations 


// //////////////////////////////////////////////////////////////////////////////
// operations with DMatrix as lhs

// gm_lhs = op(M_rhs).
template <class M>
void assign(DMatrix* gm_lhs, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
  gm_lhs->resize(	 rows(M_rhs.rows(),M_rhs.cols(),trans_rhs)
          ,cols(M_rhs.rows(),M_rhs.cols(),trans_rhs) );
  (*gm_lhs) = 0.0;
  Mp_StM(&(*gm_lhs)(),1.0,M_rhs,trans_rhs);
}

// gm_lhs = alpha * op(M_rhs).
template <class M>
void M_StM(DMatrix* gm_lhs, value_type alpha, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
  gm_lhs->resize(	 rows(M_rhs.rows(),M_rhs.cols(),trans_rhs)
          ,cols(M_rhs.rows(),M_rhs.cols(),trans_rhs) );
  (*gm_lhs) = 0.0;
  Mp_StM(&(*gm_lhs)(),alpha,M_rhs,trans_rhs);
}

// gm_lhs = op(M1_rhs1) + op(M2_rhs2).
template <class M1, class M2>
void M_MpM(DMatrix* gm_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
            ,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
  gm_lhs->resize(	 rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
          ,cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs2) );
  (*gm_lhs) = 0.0;
  DMatrixSlice gms_lhs(*gm_lhs);
  Mp_M(&gms_lhs,M1_rhs1,trans_rhs1);
  Mp_M(&gms_lhs,M2_rhs2,trans_rhs2);
}

// gm_lhs = op(M_rhs1) - op(M_rhs2).
template <class M1, class M2>
void M_MmM(DMatrix* gm_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
            ,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
  gm_lhs->resize(	 rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
          ,cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
  (*gm_lhs) = 0.0;
  DMatrixSlice gms_lhs(*gm_lhs);
  Mp_M(&gms_lhs,M1_rhs1,trans_rhs1);
  Mp_StM(&gms_lhs,-1.0,M2_rhs2,trans_rhs2);
}

// gm_lhs = alpha * op(M_rhs1) + op(gms_rhs2).
template <class M>
void M_StMpM(DMatrix* gm_lhs, value_type alpha, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MopM_assert_sizes(	 M_rhs1.rows(),M_rhs1.cols(),trans_rhs1
            ,gms_rhs2.rows(),gms_rhs2.cols(),trans_rhs2);
  assign(gm_lhs,gms_rhs2,trans_rhs2);
  Mp_StM(&(*gm_lhs)(),alpha,M_rhs1,trans_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// operations with DMatrixSlice as lhs

// gms_lhs = op(M_rhs).
template <class M>
void assign(DMatrixSlice* gms_lhs, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
  Mp_M_assert_sizes(gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
    , M_rhs.rows(), M_rhs.cols(), trans_rhs	);
  (*gms_lhs) = 0.0;
  Mp_StM(gms_lhs,1.0,M_rhs,trans_rhs);
}

// gms_lhs = alpha * op(M_rhs).
template <class M>
void M_StM(DMatrixSlice* gms_lhs, value_type alpha, const M& M_rhs, BLAS_Cpp::Transp trans_rhs) {
  Mp_M_assert_sizes(gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
    , M_rhs.rows(), M_rhs.cols(), trans_rhs	);
  (*gms_lhs) = 0.0;
  Mp_StM(gms_lhs,alpha,M_rhs,trans_rhs);
}

// gms_lhs = op(M1_rhs1) + op(M2_rhs2).
template <class M1, class M2>
void M_MpM(DMatrixSlice* gms_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
            ,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
  assert_gms_lhs(*gms_lhs, rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
               , cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
  (*gms_lhs) = 0.0;
  Mp_M(gms_lhs,M1_rhs1,trans_rhs1);
  Mp_M(gms_lhs,M2_rhs2,trans_rhs2);
}

// gms_lhs = op(M_rhs1) - op(M_rhs2).
template <class M1, class M2>
void M_MmM(DMatrixSlice* gms_lhs, const M1& M1_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MopM_assert_sizes(	 M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1
            ,M2_rhs2.rows(),M2_rhs2.cols(),trans_rhs2 );
  assert_gms_lhs(*gms_lhs, rows(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1)
               , cols(M1_rhs1.rows(),M1_rhs1.cols(),trans_rhs1) );
  (*gms_lhs) = 0.0;
  Mp_M(gms_lhs,M1_rhs1,trans_rhs1);
  Mp_StM(gms_lhs,-1.0,M2_rhs2,trans_rhs2);
}

// gms_lhs = alpha * op(M_rhs1) + op(gms_rhs2).
template <class M>
void M_StMpM(DMatrixSlice* gms_lhs, value_type alpha, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MopM_assert_sizes(	 M_rhs1.rows(),M_rhs1.cols(),trans_rhs1
            ,gms_rhs2.rows(),gms_rhs2.cols(),trans_rhs2);
  assign(gms_lhs,gms_rhs2,trans_rhs2);
  Mp_StM(gms_lhs,alpha,M_rhs1,trans_rhs1);
}

// //////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////// /////
// Level 2 BLAS

// //////////////////////////////////////////////////////////////////////////////
// += operations

// //////////////////////////////////////////////////////////////////////////////
// operations with DVector as lhs

// v_lhs = alpha * op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_StMtV(DVector* v_lhs, value_type alpha, const M& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const V& V_rhs2)
{
  MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
  v_lhs->resize(rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1));
  Vp_StMtV(&(*v_lhs)(),alpha,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// v_lhs = op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_MtV(DVector* v_lhs, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const V& V_rhs2)
{
  MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
  v_lhs->resize(rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1));
  Vp_StMtV(&(*v_lhs)(),1.0,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// //////////////////////////////////////////////////////////////////////////////
// operations with DVectorSlice as lhs

// vs_lhs = alpha * op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_StMtV(DVectorSlice* vs_lhs, value_type alpha, const M& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const V& V_rhs2)
{
  MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
  Vp_V_assert_sizes( vs_lhs->dim(), rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1) );
  Vp_StMtV(vs_lhs,alpha,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// vs_lhs = op(M_rhs1) * V_rhs2.
template <class M, class V>
void V_MtV(DVectorSlice* vs_lhs, const M& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const V& V_rhs2)
{
  MtV_assert_sizes(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1,V_rhs2.dim());
  Vp_V_assert_sizes( vs_lhs->dim(), rows(M_rhs1.rows(),M_rhs1.cols(),trans_rhs1) );
  Vp_StMtV(vs_lhs,1.0,M_rhs1,trans_rhs1,V_rhs2,0.0);
}

// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
// Level 3 BLAS

// //////////////////////////////////////////////////////////////////////////////
// += operations 

// //////////////////////////////////////////////////////////////////////////////
// = operations with DMatrix as lhs

// gm_lhs = alpha * op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_StMtM(DMatrix* gm_lhs, value_type alpha, const M1& M1_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
            , M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
  gm_lhs->resize(	  rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
          , cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
  Mp_StMtM(&(*gm_lhs)(),alpha,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0.0);
}

// gm_lhs = op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_MtM(DMatrix* gm_lhs, const M1& M1_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
            , M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
  gm_lhs->resize(	  rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
          , cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
  Mp_StMtM(&(*gm_lhs)(),1.0,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0.0);
}

// //////////////////////////////////////////////////////////////////////////////
// = operations with DMatrixSlice as lhs

// gms_lhs = alpha * op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const M1& M1_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
            , M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
  assert_gms_lhs(	  *gms_lhs
          , rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
          , cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
  Mp_StMtM(gms_lhs,alpha,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0.0);
}

// gms_lhs = op(M1_rhs1) * op(M2_rhs2).
template <class M1, class M2>
void M_MtM(DMatrixSlice* gms_lhs, const M1& M1_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  MtM_assert_sizes(	  M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
            , M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
  assert_gms_lhs(	  gms_lhs
          , rows(M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1)
          , cols(M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2) );
  Mp_StMtM(gms_lhs,1.0,M1_rhs1,trans_rhs1,M2_rhs2,trans_rhs2,0,0);
}

} // end namespace LinAlgOpPack


#endif // LIN_ALG_OP_PACK_DEF_H
