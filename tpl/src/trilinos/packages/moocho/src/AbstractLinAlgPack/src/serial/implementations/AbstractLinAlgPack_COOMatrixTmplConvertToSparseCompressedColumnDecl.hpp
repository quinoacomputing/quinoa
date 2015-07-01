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

#ifndef COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H
#define COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** @name {\bf Conversion to Fortran compatable sparse compressed column 
  * operations for COOMatrixTemplateInterface (Level 2,3 BLAS)}.
  *
  * See the ConvertToCSC class.
  */
//@{

/** \brief . */
template<class T_COOM>
size_type COOM_num_in_column(
    const T_COOM&						m
  , BLAS_Cpp::Transp					trans
  , size_type							col_offset
  , const IVector::value_type*		col_perm
  , size_type*						num_in_col	);

/** \brief . */
template<class T_COOM>
void COOM_insert_nonzeros(
    const T_COOM&						m
  , BLAS_Cpp::Transp					trans
  , value_type						alpha
  , size_type							row_offset
  , size_type							col_offset
  , const IVector::value_type*		row_perm
  , const IVector::value_type*		col_perm
  , size_type*						next_nz_in_col
  , FortranTypes::f_dbl_prec*			D_val
  , FortranTypes::f_int*				D_row_i			);

/** \brief . */
template<class T_COOM>
value_type COOM_insert_scaled_nonzeros(
    const T_COOM&						m
  , BLAS_Cpp::Transp					trans
  , value_type						scaled_max_ele
  , size_type							row_offset
  , size_type							col_offset
  , const IVector::value_type*		row_perm
  , const IVector::value_type*		col_perm
  , size_type*						next_nz_in_col
  , FortranTypes::f_dbl_prec*			D_val
  , FortranTypes::f_int*				D_row_i			);

//@}

} // end namespace AbstractLinAlgPack

#endif	// COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H
