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

#ifndef CAST_IQ_MEMBER_H
#define CAST_IQ_MEMBER_H

#include <limits.h>

#include <typeinfo>

#include "IterationPack_AlgorithmState.hpp"
#include "IterationPack_IterQuantityAccess.hpp"

namespace IterationPack {

/** \brief Base class for some of the implementation features of \c CastIQMember.
  *
  * This class is included to avoid code blot with the templates.
  */
class CastIQMemberBase {
public:
  /** \brief Name returns the name of the iteration quantity
   */
  const std::string& iq_name() const;
  /** \brief Returns if the iteration quantity exists in the state object.
   */
  bool exists_in( const AlgorithmState& s ) const;
protected:
  /** \brief . */
  CastIQMemberBase( const std::string iq_name );
  /** \brief . */
  void cache_iq_id( const AlgorithmState& s ) const;
  /** \brief . */
  void throw_cast_error( const AlgorithmState::iq_id_type iq_id, const std::string& iqa_name ) const;
  /** \brief . */
  const std::string					iq_name_;
  /** \brief . */
  mutable AlgorithmState::iq_id_type	iq_id_;
private:
  enum { NOT_SET_YET = AlgorithmState::DOES_NOT_EXIST - 1 };
  CastIQMemberBase(); // not defined and not to be called.
};	// end class CastIQMemberBase

/** \brief Template class to be used to lookup an interation quantity,
 * cast it to an <tt>\ref IterQuantityAccess "IterQuantityAccess<T>"</tt> object
 * and cache the \c iq_id for fast access later.
 *
 * The idea is that a Step class can create a data member of this type and then
 * access and interation quantity and have the \c iq_id looked up the first time.
 * 
 * The best way to use this class to access an iteration quantity is for a header
 * file to be created for the iteration quantity (or several iteration quantities)
 * that contains a new class.  For example, suppose we have an \c double object that
 * we want to add as an iteration quantity with the name "x_step".  For this we might
 * create an header file like:
 \code
    // /////////////////////////////////////////////////////////////////
    // x_step_iter_quantity.h
    
    #include "IterationPack_CastIQMember.hpp"
    
    class x_step_iq_member : public CastIQMember<double> {
    public:
        x_step_iq_member() : CastIQMember<double>("x_step") {}
    }
 \endcode 
 * Now lets suppose we have two step classes that need to access this
 * iteration quantity.  These step classes would each include
 * a data member of this new class.  For example, these
 * step classes might be implemented as:
 \code
    // /////////////////////////////////////////////////////////////////
    // MyStep1.h
    
    #include "x_step_iter_quantity.h"
    
    class MyStep1 : public AlgorithmStep {
    public:
        bool do_step( algo, ... )
        {
          AlgorithmState &s = algo.state();
          x_step_(s).set_k(0) = 5.0;
        }
    private:
      x_step_iq_member x_step_;
    }
  \endcode
  \code
    // /////////////////////////////////////////////////////////////////
    // MyStep2.h
    
    #include "x_step_iter_quantity.h"
    
    class MyStep2 : public AlgorithmStep {
    public:
        bool do_step( algo, ... )
        {
          AlgorithmState &s = algo.state();
          double x_step = x_step_(s).get_k(0);
          cout << "\nx_step = " << x_step << std::endl;
        }
    private:
      x_step_iq_member x_step_;
    }
 \endcode
 * In the above example, an <tt>O(s.num_iter_quantities())</tt> search for the \c iq_id
 * would only be performed the first time <tt>x_step_(s)</tt> was called by each
 * step object.  In later iterations, the cached \c iq_id would be used to access
 * the iteration quantity and the only price one would pay above a few
 * \a O(1) function calls in is an \a O(1) dynamic cast.
 *
 * The default constructor is not allowed by the default copy constructors and
 * assignment operators are allowed since they have the correct semantics.
 */
template < class T >
class CastIQMember : public CastIQMemberBase {
public:
  /// Construct with the name of an iteration quantity.
  CastIQMember( const std::string iq_name );
  /** \brief Get the iteration quantity from an AlgorithmState object.
    *
    * If the iteration quantity of the name iq_namt does not
    * exist then a AlgorithmState::DoesNotExist exception
    * will be thrown.  If the type of the iteration quantity
    * is not of the type IterQuantityAcess<T> (as determined
    * by dynamic_cast<T>) then the exception InvalidTypeCastException:
    * will be thrown with a helpful error message.
    */
  IterQuantityAccess<T>& operator()( AlgorithmState& s ) const;
  /** \brief . */
  const IterQuantityAccess<T>& operator()( const AlgorithmState& s ) const;
private:
  CastIQMember();	// not defined and not to be called
};	// end class CastIQMember<T>

// //////////////////////////////////////////
// Definition of template members

template < class T >
CastIQMember<T>::CastIQMember( const std::string iq_name )
  :  CastIQMemberBase(iq_name)
{}

template < class T >
IterQuantityAccess<T>&
CastIQMember<T>::operator()( AlgorithmState& s ) const
{
  cache_iq_id(s);
  if( iq_id_ == AlgorithmState::DOES_NOT_EXIST )
    throw_cast_error(iq_id_,TypeNameTraits<T>::name());
  IterQuantityAccess<T>
    *p = dynamic_cast<IterQuantityAccess<T>*>( &s.iter_quant( iq_id_ ) );
  if( !p )
    throw_cast_error(iq_id_,TypeNameTraits<T>::name());
  return *p;	
}

template < class T >
const IterQuantityAccess<T>&
CastIQMember<T>::operator()( const AlgorithmState& s ) const
{
  cache_iq_id(s);
  if( iq_id_ == AlgorithmState::DOES_NOT_EXIST )
    throw_cast_error(iq_id_,TypeNameTraits<T>::name());
  const IterQuantityAccess<T>
    *p = dynamic_cast<const IterQuantityAccess<T>*>( &s.iter_quant( iq_id_ ) );
  if( !p )
    throw_cast_error(iq_id_,TypeNameTraits<T>::name());
  return *p;	
}

}	// namespace IterationPack

#endif	// CAST_IQ_MEMBER_H
