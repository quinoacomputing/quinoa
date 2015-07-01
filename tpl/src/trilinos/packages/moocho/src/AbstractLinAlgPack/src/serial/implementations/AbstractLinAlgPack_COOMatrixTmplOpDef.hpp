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

#ifndef COO_MATRIX_TMPL_OP_DEF_H
#define COO_MATRIX_TMPL_OP_DEF_H

#include "AbstractLinAlgPack_COOMatrixTmplOpDecl.hpp"

#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace AbstractLinAlgPack {

using BLAS_Cpp::trans_not;

using DenseLinAlgPack::Mp_M_assert_sizes;
using DenseLinAlgPack::Vp_MtV_assert_sizes;
using DenseLinAlgPack::Mp_MtM_assert_sizes;

// gms_lhs += alpha * coom_rhs (time = O(coom_rhs.nz()), space = O(1))
template<class T_COOM>
void Mp_StCOOM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs
  , BLAS_Cpp::Transp trans_rhs)
{
  Mp_M_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
            , coom_rhs.rows(), coom_rhs.cols(), trans_rhs			);
  typename T_COOM::difference_type
    i_o	= coom_rhs.row_offset(),
    j_o	= coom_rhs.col_offset();
  if(trans_rhs == BLAS_Cpp::no_trans)
    for(typename T_COOM::const_iterator itr = coom_rhs.begin(); itr != coom_rhs.end(); ++itr)
      (*gms_lhs)(itr->row_i()+i_o,itr->col_j()+j_o) += alpha * itr->value();
  else
    for(typename T_COOM::const_iterator itr = coom_rhs.begin(); itr != coom_rhs.end(); ++itr)
      (*gms_lhs)(itr->col_j()+j_o,itr->row_i()+i_o) += alpha * itr->value();
}

// vs_lhs += alpha * coom_rhs1 * vs_rhs2 (BLAS xGEMV) (time = O(coom_rhs.nz()), space = O(1))
template<class T_COOM>
void Vp_StCOOMtV(DVectorSlice* vs_lhs, value_type alpha, const T_COOM& coom_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2)
{
  Vp_MtV_assert_sizes( vs_lhs->dim(), coom_rhs1.rows(), coom_rhs1.cols(), trans_rhs1, vs_rhs2.dim() );
  typename T_COOM::difference_type
    i_o	= coom_rhs1.row_offset(),
    j_o	= coom_rhs1.col_offset();
  if(trans_rhs1 == BLAS_Cpp::no_trans)
    for(typename T_COOM::const_iterator itr = coom_rhs1.begin(); itr != coom_rhs1.end(); ++itr)
      (*vs_lhs)(itr->row_i()+i_o) += alpha * itr->value() * vs_rhs2(itr->col_j()+j_o);
  else
    for(typename T_COOM::const_iterator itr = coom_rhs1.begin(); itr != coom_rhs1.end(); ++itr)
      (*vs_lhs)(itr->col_j()+j_o) += alpha * itr->value() * vs_rhs2(itr->row_i()+i_o);
}

namespace UtilityPack {

// op(gms_lhs) += alpha * op(gms_rhs1) * op(coom_rhs2) (BLAS xGEMM)
template<class T_COOM>
void imp_Mp_StMtCOOM(DMatrixSlice& gms_lhs, BLAS_Cpp::Transp trans_lhs, value_type alpha
  , const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2 );

}	// end namespace UtilityPack


// gms_lhs += alpha * op(coom_rhs1) * op(gms_rhs2) (right) (BLAS xGEMM)
template<class T_COOM>
void Mp_StCOOMtM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_MtM_assert_sizes(  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
            , coom_rhs1.rows(), coom_rhs1.cols(), trans_rhs1
            , gms_rhs2.rows() ,gms_rhs2.cols(), trans_rhs2			);
  UtilityPack::imp_Mp_StMtCOOM( gms_lhs, BLAS_Cpp::trans, alpha, gms_rhs2
    , trans_not(trans_rhs2), coom_rhs1, trans_not(trans_rhs1)				);
}

// gms_lhs += alpha * op(gms_rhs1) * op(coom_rhs2) (left) (BLAS xGEMM)
template<class T_COOM>
void Mp_StMtCOOM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_MtM_assert_sizes(  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
            , gms_rhs1.rows() ,gms_rhs1.cols(), trans_rhs1
            , coom_rhs2.rows(), coom_rhs2.cols(), trans_rhs2		);
  UtilityPack::imp_Mp_StMtCOOM( gms_lhs, BLAS_Cpp::no_trans, alpha, gms_rhs1
    , trans_rhs1, coom_rhs2, trans_rhs2										);
}

// ToDo: implement the level-3 BLAS operations below

// gms_lhs = alpha * op(coom_rhs1) * op(sym_rhs2) (right) (BLAS xSYMM)
//template<class T_COOM>
//void Mp_StCOOMtSM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceSym& sym_rhs2, BLAS_Cpp::Transp trans_rhs2);

// gms_lhs = alpha * op(sym_rhs1) * op(coom_rhs2) (left) (BLAS xSYMM)
//template<class T_COOM>
//void Mp_StSMtCOOM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceSym& sym_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2);

// gms_lhs = alpha * op(coom_rhs1) * op(tri_rhs2) (right) (BLAS xTRMM)
//template<class T_COOM>
//void Mp_StCOOMtSM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceTri& tri_rhs2, BLAS_Cpp::Transp trans_rhs2);

// gms_lhs = alpha * op(tri_rhs1) * op(coom_rhs2) (left) (BLAS xTRMM)
//template<class T_COOM>
//void Mp_StSMtCOOM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceTri& tri_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2);

namespace UtilityPack {

// op(gms_lhs) += alpha * op(gms_rhs1) * op(coom_rhs2) (BLAS xGEMM)
//
// This function perform this operation by looping through the nonzero elements of
// coom_rhs2 and for each nonzero element (val,i,j) it performs the following operation.
//
//	op(gms_lhs).col(j) += (alpha * val) * op(gms_rhs1).col(i)
//
//
//			jth								ith						jth
//		[	#			]		[			#		]			[				]
//		[	#			]		[			#		]			[				]
//	m	[	#			]	+=	[			#		]	*		[				] n
//		[	#			]		[			#		]		ith	[	val			]
//				n						p						[				]
//																		p
//			op(gms_lhs)				op(gms_rhs1)				   op(coom_rhs2)
//
// The number of arithmetic operations performed are:	
//	floats = (2*m + 1) * nz
//
// Strictly speeking the number of memory references is:
//	mem_refs = (2*m + 1) * nz
// but this does not take into account that elements are accessed by columns
// and this has some ramifications on cache effects and paging.  If op(gms_lhs)
// == gms_lhs' or op(gms_rhs1) == gms_rhs1' then elements in a column are
// not adjacent to each other and if m is large enough each element may
// even reside on a seperate page of memory.  On Win32 0x86 systems a page is
// 4 K so 4,000 (bytes/page) / 8 (bytes/double) = 500 doubles / page.  If
// the working set of pages is small this could cause some serious page thrashing
// for large m.
//
// Another concideration is to sorting order of the elements in the COO matrix.
// If op(coom_rhs2) is sorted by row then columns of op(gms_lhs) will be accessed
// consecutivly and will result in better performance.  The same goes for op(gms_rhs1)
// if op(coom_rhs2) is sorted by column.
//
// There is opertunity for some vectorization and it is handled by calling
// DenseLinAlgPack::Vp_StV(...).
//
template<class T_COOM>
void imp_Mp_StMtCOOM(DMatrixSlice* gms_lhs, BLAS_Cpp::Transp trans_lhs, value_type alpha
  , const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  using DenseLinAlgPack::col;

  typename T_COOM::difference_type
    i_o	= coom_rhs2.row_offset(),
    j_o	= coom_rhs2.col_offset();
  for(typename T_COOM::const_iterator itr = coom_rhs2.begin(); itr != coom_rhs2.end(); ++itr) {
    size_type	i	= rows( itr->row_i() + i_o , itr->col_j() + j_o , trans_rhs2 ),
          j	= cols( itr->row_i() + i_o , itr->col_j() + j_o , trans_rhs2 );
    //	op(gms_lhs).col(j) += (alpha * val) * op(gms_rhs1).col(i)
    DenseLinAlgPack::Vp_StV(	&col(*gms_lhs,trans_lhs,j), alpha * itr->value()
      , col(gms_rhs1,trans_rhs1,i) );
  }
}

}	// end namespace UtilityPack

} // end namespace AbstractLinAlgPack

#endif	// COO_MATRIX_TMPL_OP_DEF_H
