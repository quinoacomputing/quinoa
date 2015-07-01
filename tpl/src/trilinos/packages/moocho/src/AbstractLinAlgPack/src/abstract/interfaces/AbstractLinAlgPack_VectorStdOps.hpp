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

#ifndef ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H
#define ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H

#include "AbstractLinAlgPack_VectorMutable.hpp"

namespace AbstractLinAlgPack {

/** \defgroup VectorStdOps_grp Collection of standard vector operations.
 */
//@{

/** \defgroup VectorStdOps_ROp_grp Reduction operations */
//@{

/** \brief result = sum( v_rhs(i), i = 1,,,dim )
 */
value_type sum( const Vector& v_rhs );

/** \brief result = v_rhs1' * v_rhs2
 */
value_type dot( const Vector& v_rhs1, const Vector& v_rhs2 );

/** \brief result = v_rhs1' * sv_rhs2
 */
value_type dot( const Vector& v_rhs1, const SpVectorSlice& sv_rhs2 );

/** \brief result = sv_rhs1' * v_rhs2
 */
value_type dot( const SpVectorSlice& sv_rhs1, const Vector& v_rhs2 );

/** \brief Compute the maximum element in a vector.
 *
 * @param  v        [in] The vector being searched
 * @param  max_v_j  [out] The value of the element with the max abs value.
 * @param  max_j    [out] The index of the element with the max abs value.
 *
 * Returns:
 \verbatim

 max_v_j = v(max_j) s.t. |v(max_j)| <= |v(j)|, for j = 1...n
 \endverbatim
 * If there is a tie, the lowest index is returned so that the
 * result is unique no matter what order the vector elements are
 * searched.
 */
void max_abs_ele( const Vector& v, value_type* max_v_j, index_type* max_j ); 

//@}

/** \defgroup VectorStdOps_TOp_grp Transformation operations */
//@{

/** \brief v_lhs += alpha
 */
void Vp_S( VectorMutable* v_lhs, const value_type& alpha );

/** \brief v_lhs *= alpha
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>v_lhs = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 */
void Vt_S( VectorMutable* v_lhs, const value_type& alpha );

/** \brief v_lhs = alpha * v_rhs + v_lhs
 */
void Vp_StV( VectorMutable* v_lhs, const value_type& alpha, const Vector& v_rhs );

/** \brief v_lhs = alpha * sv_rhs + v_lhs
 */
void Vp_StV( VectorMutable* v_lhs, const value_type& alpha, const SpVectorSlice& sv_rhs );

/** \brief v_lhs(i) += alpha * v_rhs1(i) * v_rhs2(i), i = 1,,,dim.
 */
void ele_wise_prod(
  const value_type& alpha, const Vector& v_rhs1, const Vector& v_rhs2
  ,VectorMutable* v_lhs );

/** \brief v_lhs(i) = alpha * v_rhs1(i) / v_rhs2(i), i = 1,,,dim.
 */
void ele_wise_divide(
  const value_type& alpha, const Vector& v_rhs1, const Vector& v_rhs2
  ,VectorMutable* v_lhs );

/** \brief Seed the random number generator
  */
void seed_random_vector_generator( unsigned int );

/** \brief Generate a random vector with elements uniformly
  * distrubuted elements.
  * 
  * The elements are randomly generated between <tt>[l,u]</tt>.
  */
void random_vector( value_type l, value_type u, VectorMutable* v );

/** \brief Compute the sign of each element in an input vector.
 *
 \verbatim

         / -1.0 : if v(i)  < 0.0
 z(i) =  |  0.0 : if v(i) == 0.0
         \ +1.0 : if v(i)  < 0.0

 , for i = 1...n
 \endverbatim
 */
void sign(
  const Vector      &v
  ,VectorMutable    *z
  );

//@}

//@}

} // end namespace AbstractLinAlgPack

// /////////////////////////////////////
// Inline members

inline
AbstractLinAlgPack::value_type
AbstractLinAlgPack::dot( const SpVectorSlice& sv_rhs1, const Vector& v_rhs2 )
{
  return dot(v_rhs2,sv_rhs1);
}

#endif // ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H
