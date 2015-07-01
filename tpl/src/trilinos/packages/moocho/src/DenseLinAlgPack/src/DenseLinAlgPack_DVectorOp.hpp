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

// See DenseLinAlgPack_DMatrixOp.hpp for description of naming convensions

#ifndef VECTOROP_H
#define VECTOROP_H

#include "DenseLinAlgPack_DVectorAssign.hpp"

/* * @name {\bf Basic DVector Operation Functions (Level-1 BLAS)}.
  *
  * These are functions that perform basic operations with vectors such as element-wise
  * linear algebra operations (e.g. v1 = v2 + v3, v1 = sin(v2)) and other vector 
  * related functions.  The functions that have vectors as lhs arguments come in
  * two varieties: those with a DVector object as the lhs argument (v_lhs), those with
  * a DVectorSlice object as a lhs argument (vs_lhs).  Having different functions
  * for DVector and DVectorSlice objects as lhs arguments is important because
  * Vectors resize to the rhs expression while DVectorSlice objects do not.
  *
  * Only DVectorSlice objects are used as rhs arguments however.  When DVector objects
  * are used to call these fucntions as rhs arguments, the implicit type conversion
  * to a const temp DVectorSlice will be performed to make the call work.
  * 
  * The implementations of these functions takes care of the following details:
  *
  * <ul>
  *	<li> Resizing DVector LHS on assignment
  *	<li> Test for aliasing of assign(...) but not other functions
  *	<li> Check preconditions (sizes of arguments) if LINALGPACK_CHECK_RHS_SIZES is defined
  * </ul>
  *
  * These functions share common behavior and precondtions which are listed below.
  *
  * Preconditions for functions with a DVectorSlice object (vs_lhs) as a lhs argument
  * (e.g. vs_lhs = abs(vs_rhs), vs_lhs = vs_rhs1 + vs_rhs2).
  *	<ul>
  * <li> #vs_lhs.size() ==# size of rhs expression  (throw #std::length_error#)
  * </ul>
  *
  * Preconditions for functions with two DVectorSlice objects (vs_rhs1, vs_rhs2) rhs arguments
  * (e.g. v_lhs = pow(vs_rhs1,vs_rhs2), result = trans(vs_rhs1,vs_rhs2)):
  *	<ul>
  * <li> #vs_rhs1.size() == vs_rhs2.size()#  (throw #std::length_error#)
  * </ul>
  *
  * Algebric functions are named according to the types of their arguments.  For example,
  * the function for the operation vs_lhs = vs_rhs1 - vs_rhs2 is named V_VmV(...).
  * For a description of this namming format see
  * \Ref{LinAlgOpPack}
  */

// @{
//		begin Basic DVector Operation Functions

namespace DenseLinAlgPack {

/* * @name {\bf Algebraic Functions}.
  *
  * The functions assign(...) are used by the implementation of the assignment operators for
  * DVector and DVectorSlice and therefore the user can use the assignment operator to
  * perform the copies.
  */

// @{
//		begin Algebraic Functions


/// vs_lhs += alpha
void Vp_S(DVectorSlice* vs_lhs, value_type alpha);
/// vs_lhs *= alpha (BLAS xSCAL) (*** Note that alpha == 0.0 is handeled as vs_lhs = 0.0)
void Vt_S(DVectorSlice* vs_lhs, value_type alpha);
/// vs_lhs += alpha * vs_rhs (BLAS xAXPY)
void Vp_StV(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs);

/// v_lhs = alpha (elementwise)
//void assign(DVector* v_lhs, value_type alpha);
/// v_lhs = vs_rhs.
//void assign(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = vs_rhs1 + vs_rhs2
void V_VpV(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// v_lhs = vs_rhs1 - vs_rhs2
void V_VmV(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// v_lhs = - vs_rhs
void V_mV(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = alpha * vs_rhs
void V_StV(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs);

/// vs_lhs = alpha (elementwise)
//void assign(DVectorSlice* vs_lhs, value_type alpha);
/// vs_lhs = vs_rhs
//void assign(DVectorSlice* vs_lhs, const DVectorSlice& vx_rhs);
/// vs_lhs = vs_rhs1 + vs_rhs2
void V_VpV(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// vs_lhs = vs_rhs1 - vs_rhs2
void V_VmV(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// vs_lhs = - vs_rhs
void V_mV(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = alpha * vs_rhs
void V_StV(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs);

/** \brief . */
/* Apply a plane (Givens) rotation.
 *
 * [  c  s ] * [ x' ] -> [ x' ]
 * [ -s  c ]   [ y' ]    [ y' ]
 *
 * See "Handbook for Matrix Computations" section 2.4
 */
void rot( const value_type c, const value_type s, DVectorSlice* x, DVectorSlice* y );

//		end Algebraic Functions
// @}

/* * @name {\bf Elementwise Math DVector / DVectorSlice Functions}. */

// @{
//		begin Elementsize Math Functions

/// vs_lhs = abs(vs_rhs)
void abs(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = asin(vs_rhs)
void asin(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = acos(vs_rhs)
void acos(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = atan(vs_rhs)
void atan(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = atan(vs_rhs1/vs_rhs2)
void atan2(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// vs_lhs = atan(vs_rhs/alpha)
void atan2(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs, value_type alpha);
/// vs_lhs = atan(alpha/vs_rhs)
void atan2(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs);
/// vs_lhs = cos(vs_rhs)
void cos(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = cosh(vs_rhs)
void cosh(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = exp(vs_rhs)
void exp(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = max(vs_rhs1,vs_rhs2)
void max(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// vs_lhs = max(alpha,vs_rhs)
void max(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs);
/// vs_lhs = min(vs_rhs1,vs_rhs2)
void min(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// vs_lhs = mim(alpha,vs_rhs)
void min(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs);
/// vs_lhs = pow(vs_rhs1,vs_rhs2)
void pow(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// vs_lhs = pow(vs_rhs,alpha)
void pow(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs, value_type alpha);
/// vs_lhs = pow(vs_rhs,n) 
void pow(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs, int n);
/// vs_lhs = pow(alpha,vs_rhs)
void pow(DVectorSlice* vs_lhs, value_type alpha, const DVectorSlice& vs_rhs);
/// vs_lhs(i) = vs_rhs1(i) * vs_rhs2(i), i = 1...n
void prod(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// vs_lhs = sqrt(vs_rhs)
void sqrt(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = sin(vs_rhs)
void sin(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = sinh(vs_rhs)
void sinh(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = tan(vs_rhs)
void tan(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = tanh(vs_rhs)
void tanh(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);

/// v_lhs = abs(vs_rhs)
void abs(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = asin(vs_rhs)
void asin(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = acos(vs_rhs)
void acos(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = atan(vs_rhs)
void atan(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = atan(vs_rhs1/vs_rhs2)
void atan2(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// v_lhs = atan(vs_rhs/alpha)
void atan2(DVector* v_lhs, const DVectorSlice& vs_rhs, value_type alpha);
/// v_lhs = atan(alpha/vs_rhs)
void atan2(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs);
/// v_lhs = cos(vs_rhs)
void cos(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = cosh(vs_rhs)
void cosh(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = exp(vs_rhs)
void exp(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = max(vs_rhs1,vs_rhs2)
void max(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// v_lhs = max(alpha,vs_rhs)
void max(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs);
/// v_lhs = min(vs_rhs1,vs_rhs2)
void min(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// v_lhs = mim(alpha,vs_rhs)
void min(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs);
/// v_lhs = pow(vs_rhs1,vs_rhs2)
void pow(DVector* v_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// v_lhs = pow(vs_rhs,alpha)
void pow(DVector* v_lhs, const DVectorSlice& vs_rhs, value_type alpha);
/// v_lhs = pow(vs_rhs,n) 
void pow(DVector* v_lhs, const DVectorSlice& vs_rhs, int n);
/// v_lhs = pow(alpha,vs_rhs)
void pow(DVector* v_lhs, value_type alpha, const DVectorSlice& vs_rhs2);
/// v_lhs = sqrt(vs_rhs)
void sqrt(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = sin(vs_rhs)
void sin(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = sinh(vs_rhs)
void sinh(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = tan(vs_rhs)
void tan(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs = tanh(vs_rhs)
void tanh(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// v_lhs(i) = vs_rhs1(i) * vs_rhs2(i), i = 1...n
void prod( DVector* vs_lhs, const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2 );

//		end Elementsize Math Functions
// @}

/* * @name {\bf Scalar Returning and Misc DVectorSlice Functions}. */

// @{
//		begin Scalar Returning DVectorSlice Functions}

/// result = vs_rhs1' * vs_rhs2 (BLAS xDOT)
value_type dot(const DVectorSlice& vs_rhs1, const DVectorSlice& vs_rhs2);
/// result = max(vs_rhs)
value_type max(const DVectorSlice& vs_rhs);
/// result = min(vs_rhs)
value_type min(const DVectorSlice& vs_rhs);
/// result = ||vs_rhs||1 (BLAS xASUM)
value_type norm_1(const DVectorSlice& vs_rhs);
/// result = ||vs_rhs||2 (BLAS xNRM2)
value_type norm_2(const DVectorSlice& vs_rhs);
/// result = ||vs_rhs||infinity (BLAS IxAMAX)
value_type norm_inf(const DVectorSlice& vs_rhs);

// Misc. operations

/// swap(vs1, vs2). Swaps the contents of vs1 and vs2
void swap(DVectorSlice* vs1, DVectorSlice* vs2);

//		end Scalar Returning DVectorSlice Functions
// @}

} // end namespace DenseLinAlgPack

//		end Basic DVector Operation Functions
// @}

#endif // VECTOROP_H
