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

#ifndef GIP_CAST_IQ_H
#define GIP_CAST_IQ_H

#include <stdexcept>
#include <typeinfo>

#include "IterationPack_AlgorithmState.hpp"
#include "IterationPack_IterQuantityAccess.hpp"

namespace IterationPack {

/** \brief Lookup an iteration quantity by name and cast it
  * to an <tt>IterQuantityAccess<T></tt> of the given type \c T.
  * If the iteration quantity of that name does not
  * exist then a <tt>AlgorithmState::DoesNotExist</tt> exception
  * will be thrown.  If the type of the iteration quantity
  * is not of the type <tt>IterQuantityAcess<T></tt> (as determined
  * by <tt>dynamic_cast<T></tt>) then the exception \c InvalidTypeCastException:
  * will be thrown with a helpful error message.
  *
  * Note that using this function always cost <tt>O(s.num_iter_quantities())</tt>
  * everytime it is called.  Therefore, the developer should consider using the
  * class \c CastIQMember instead if it is appropriate.
  */
template<class T>
IterQuantityAccess<T>& cast_iq(
  AlgorithmState& state, const std::string& iq_name );

/** \brief . */
template<class T>
const IterQuantityAccess<T>& cast_iq(
  const AlgorithmState& state, const std::string& iq_name );

/** \brief Lookup an iteration quantity using its id and cast it
 * to an <tt>IterQuantityAccess<T></tt> of the given type \c T.
 *
 * @param  state  [in] 
 *
 * If the iteration quantity with that id does not
 * exist then a <tt>AlgorithmState::DoesNotExist</tt> exception
 * will be thrown.  If the type of the iteration quantity
 * is not of the type <tt>IterQuantityAcess<T></tt> (as determined
 * by <tt>dynamic_cast<T></tt>) then the exception \c InvalidTypeCastException:
 * will be thrown with a helpful error message.
 *
 * Note that using this function always cost just <tt>O(1)</tt> everytime
 * it is called.
 */
template<class T>
IterQuantityAccess<T>& cast_iq(
  AlgorithmState& state, const AlgorithmState::iq_id_type iq_id, const std::string& iq_name );

/** \brief . */
template<class T>
const IterQuantityAccess<T>& cast_iq(
  const AlgorithmState& state, const AlgorithmState::iq_id_type iq_id, const std::string& iq_name );

// Helper function

void imp_cast_iq_throw_error(
  const std::string&                 iq_name
  ,const std::string&                iq_is_type_name
  ,const std::string&                iq_want_type_name
  );

void imp_cast_iq_throw_error(
  const AlgorithmState::iq_id_type   iq_id
  ,const std::string&                iq_name
  ,const std::string&                iq_is_type_name
  ,const std::string&                iq_want_type_name
  );

// ///////////////////
// Inline definitions

template<class T>
//inline
IterQuantityAccess<T>& cast_iq(
  AlgorithmState& state, const std::string& iq_name )
{
  IterQuantity
     &iq = state.iter_quant( iq_name );
  IterQuantityAccess<T>
    *p = dynamic_cast<IterQuantityAccess<T>*>( &iq );
      // will throw exception if iq_name does not exist
  if( !p )
    imp_cast_iq_throw_error( iq_name, typeName(iq), TypeNameTraits<T>::name() );
  return *p;	
}

template<class T>
//inline
const IterQuantityAccess<T>& cast_iq(
  const AlgorithmState& state, const std::string& iq_name )
{
  const IterQuantity
     &iq = state.iter_quant( iq_name );
  const IterQuantityAccess<T>
    *p = dynamic_cast<const IterQuantityAccess<T>*>( &iq );
      // will throw exception if iq_name does not exist
  if( !p )
    imp_cast_iq_throw_error( iq_name, typeName(iq), TypeNameTraits<T>::name() );
  return *p;	
}

template<class T>
//inline
IterQuantityAccess<T>& cast_iq(
  AlgorithmState& state, const AlgorithmState::iq_id_type iq_id, const std::string& iq_name )
{
  IterQuantity
     &iq = state.iter_quant( iq_id );
  IterQuantityAccess<T>
    *p = dynamic_cast<IterQuantityAccess<T>*>( &iq );
      // will throw exception if iq_name does not exist
  if( !p )
    imp_cast_iq_throw_error( iq_id, iq_name, typeName(iq), TypeNameTraits<T>::name() );
  return *p;	
}

template<class T>
//inline
const IterQuantityAccess<T>& cast_iq(
  const AlgorithmState& state, const AlgorithmState::iq_id_type iq_id, const std::string& iq_name )
{
  const IterQuantity
     &iq = state.iter_quant( iq_id );
  const IterQuantityAccess<T>
    *p = dynamic_cast<const IterQuantityAccess<T>*>( &iq );
      // will throw exception if iq_name does not exist
  if( !p )
    imp_cast_iq_throw_error( iq_id, iq_name, typeName(iq), TypeNameTraits<T>::name() );
  return *p;
}

}	// namespace IterationPack

#endif	// GIP_CAST_IQ_H
