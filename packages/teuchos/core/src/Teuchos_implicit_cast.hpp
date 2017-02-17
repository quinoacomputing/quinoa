// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_IMPLICIT_CAST_HPP
#define TEUCHOS_IMPLICIT_CAST_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

/** \brief Perform an implicit cast of concrete types with the casted object
 * returned by value.
 *
 * This function is used as:

 \code
    TypeTo myCast( const TypeFrom& a )
    {
      return Teuchos::implicit_cast<TypeTo>(a);
    }
 \endcode

 *
 * This function will only compile for types where an implicit conversion from
 * objects of type <tt>TypeFrom</tt> to type <tt>TypeTo</tt> exists.  Note
 * that this is a weaker operation than a <tt>static_cast<TypeTo>(t)</tt> in
 * that the static cast will actually compile in cases where the implicit
 * conversion would not compile and may result in incorrect code.  This
 * function can not result in incorrect code, assuming that the implicit
 * conversions themselves do no result in incorrect code (which is another
 * matter all together).
 *
 * This function is especially helpful when one needs to be careful of what
 * specific type is passed in as a formal argument to a function and in
 * comparing values.  In particular, using this function is far safer than
 * using <tt>TypeTo(t)</tt> in cases where <tt>TypeTo</tt> is a built in type
 * since <tt>TypeTo(t)</tt> in these cases is equivalent to <tt>(TypeTo)t</tt>
 * which is an unchecked sledge-hammer cast of the worst kind.
 *
 * \ingroup teuchos_language_support_grp
 */
template<class TypeTo, class TypeFrom>
inline TypeTo implicit_cast( const TypeFrom& t ) { return t; }

/** \brief Perform an implicit cast of reference types with a reference being
 * returned.
 *
 * This function is used as:

 \code
    TypeTo& myPtrCast( TypeFrom &ref1 )
    {
      return Teuchos::implicit_ref_cast<TypeTo>(ref1);
    }
 \endcode

 * This function will only compile for types where an implicit conversion from
 * references of type <tt>TypeFrom&</tt> to references of <tt>TypeTo&</tt>
 * exists.  It is allowed for the type <tt>TypeFrom</tt> and <tt>TypeTo</tt>
 * to actually be <tt>const</tt> types.  For example, we can have <tt>TypeFrom
 * = const std::iostream</tt> and <tt>TypeTo = const std::ostream</tt> and
 * <tt>Teuchos::implicit_ref_cast<const
 * std::ostream>(Teuchos::getConst(std::cout))</tt> would compile just fine.

 * Note that this is a weaker operation than a
 * <tt>static_cast<TypeTo&>(t)</tt> in that the static cast will actually
 * compile in cases where the implicit conversion would not compile and may
 * result in incorrect code.  For example, a static cast from a base to a
 * derived class will compile (and may be wrong) while this implicit cast
 * function will not compile for casts from base to derived classes.  This
 * function can not result in incorrect code, assuming that the implicit
 * conversions themselves do no result in incorrect code (which is another
 * matter all together).
 *
 * \ingroup teuchos_language_support_grp
 */
template<class TypeTo, class TypeFrom>
inline TypeTo& implicit_ref_cast( TypeFrom& t ) { return t; }

/** \brief Perform an implicit cast of pointer types with a pointer being
 * returned.
 *
 * This function is used as:

 \code
    TypeTo* myPtrCast( TypeFrom *ptr1 )
    {
      return Teuchos::implicit_ptr_cast<TypeTo>(ptr1);
    }
 \endcode

 * This function will only compile for types where an implicit conversion from
 * pointers of type <tt>TypeFrom*</tt> to pointers of type <tt>TypeTo*</tt>
 * exists.  It is allowed for the type <tt>TypeFrom</tt> and <tt>TypeTo</tt>
 * to actually be <tt>const</tt> types.  For example, we can have <tt>TypeFrom
 * = const std::iostream</tt> and <tt>TypeTo = const std::ostream</tt> and
 * <tt>Teuchos::implicit_ptr_cast<const
 * std::ostream>(&Teuchos::getConst(std::cout))</tt> would compile just fine.
 *
 * Note that this is a weaker operation than a
 * <tt>static_cast<TypeTo*>(t)</tt> in that the static cast will actually
 * compile in cases where the implicit conversion would not compile and may
 * result in incorrect code.  For example, a static cast up from a base class
 * to a derived class will compile (and may be wrong) while this implicit cast
 * function will not compile for such casts.  This function can not result in
 * incorrect code, assuming that the implicit conversions themselves do no
 * result in incorrect code (which is another matter all together).
 *
 * \ingroup teuchos_language_support_grp
 */
template<class TypeTo, class TypeFrom>
inline TypeTo* implicit_ptr_cast( TypeFrom* t ) { return t; }

} // end namespace Teuchos

#endif // TEUCHOS_IMPLICIT_CAST_HPP
