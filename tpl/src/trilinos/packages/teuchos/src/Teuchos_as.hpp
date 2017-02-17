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

#ifndef TEUCHOS_AS_HPP
#define TEUCHOS_AS_HPP

#include "Teuchos_Assert.hpp"

#ifdef HAVE_TEUCHOS_QD
#include <qd/qd_real.h>
#include <qd/dd_real.h>
#endif

namespace Teuchos {


/** \brief Default traits class for all conversions of value types.
 *
 * This class should never be called directly by clients.  Instead, use the
 * <tt>as()</tt> and <tt>asSafe()</tt> template functions.
 *
 * This default traits class simply does an implicit type conversion.
 * Therefore, any conversions that are built into the language and are safe do
 * not need a traits class specialization and should not generate any compiler
 * warnings.  For example, the conversions <tt>float</tt> to <tt>double</tt>,
 * <tt>short type</tt> to <tt>type</tt>, <tt>type</tt> to <tt>long type</tt>,
 * and an enum value to <tt>int</tt> are all always value preserving and
 * should never result in a compiler warning or any aberrant runtime behavior.
 *
 * All other conversions that cause compiler warnings and/or could result in
 * aberrant runtime behavior (e.g. <tt>type</tt> to and from <tt>unsigned
 * type</tt>, to and from floating point and integral types, etc.), or do not
 * have compiler defined conversions (e.g. <tt>std::string</tt> to
 * <tt>int</tt>, <tt>double</tt> etc.) should be given specializations of this
 * class template.  If an unsafe or non-supported conversion is requested by
 * a client (i.e. through <tt>as()</tt> or <tt>asSafe()</tt>) then this
 * default traits class will be instantiated and the compiler will either
 * generate a warning message (if the conversion is supported but is unsafe)
 * or will not compile the code (if the conversion is not supported by default
 * in C++).  When this happens, a specialization can be added or the client
 * code can be changed to avoid the conversion.
 *
 * \ingroup teuchos_language_support_grp
 */
template<class TypeTo, class TypeFrom>
class ValueTypeConversionTraits {
public:
  static TypeTo convert( const TypeFrom t )
    {
      return t;
      // This default implementation is just an implicit conversion and will
      // generate compiler warning on dangerous conversions.
    }
  static TypeTo safeConvert( const TypeFrom t )
    {
      return t;
      // This default implementation is just an implicit conversion and will
      // generate compiler warning on dangerous conversions.  No runtime
      // checking can be done by default; only specializations can define
      // meaningful and portable runtime checks of conversions.
    }
};

/** \brief Perform an debug-enabled checked conversion from one value type
 * object to another.
 *
 * This function is used as:

 \code
    TypeTo myConversion( const TypeFrom& a )
    {
      return Teuchos::as<TypeTo>(a);
    }
 \endcode 

 * This is just an interface function what calls the traits class
 * <tt>ValueTypeConversionTraits</tt> to perform the actual conversion.  All
 * specializations of behavior is done through specializations of the
 * <tt>ValueTypeConversionTraits</tt> class (which should be done in the
 * <tt>Teuchos</tt> namespace).
 *
 * When debug checking is turned on (e.g. when the <tt>TEUCHOS_DEBUG</tt>
 * macro is defined by the <tt>--enable-teuchos-debug</tt> configure option),
 * then the checked conversion function
 * <tt>ValueTypeConversionTraits<TypeTo,TypeFrom>::safeConvert(t)</tt> is
 * called.  When debug checking is not turned on, the unchecked
 * <tt>ValueTypeConversionTraits<TypeTo,TypeFrom>::convert(t)</tt> function is
 * called.
 *
 * For cases where the checking should always be done (i.e. to validate user
 * data), use the <tt>asSafe()</tt> version of this function.
 *
 * \ingroup teuchos_language_support_grp
 */
template<class TypeTo, class TypeFrom>
inline TypeTo as( const TypeFrom& t )
{
#ifdef TEUCHOS_DEBUG
  return ValueTypeConversionTraits<TypeTo,TypeFrom>::safeConvert(t);
#else
  return ValueTypeConversionTraits<TypeTo,TypeFrom>::convert(t);
#endif
}


/** \brief Perform an always checked conversion from one value type object to
 * another.
 *
 * This function is used as:

 \code
    TypeTo mySafeConversion( const TypeFrom& a )
    {
      return Teuchos::asSafe<TypeTo>(a);
    }
 \endcode 

 * This is just an interface function what calls the traits class
 * <tt>ValueTypeConversionTraits</tt> to perform the actual conversion.  All
 * specializations of behavior is done through specializations of the
 * <tt>ValueTypeConversionTraits</tt> class (which should be done in the
 * <tt>Teuchos</tt> namespace).
 *
 * This function always calls
 * <tt>ValueTypeConversionTraits<TypeTo,TypeFrom>::safeConvert(t)</tt>
 * independent of whether <tt>TEUCHOS_DEBUG</tt> is defined or not, which
 * ensures that the conversions are always runtime checked and therefore well
 * defined.
 *
 * For cases where the checking should only be done in a debug build, use the
 * the <tt>as()</tt> version of this function.
 *
 * \ingroup teuchos_language_support_grp
 */
template<class TypeTo, class TypeFrom>
inline TypeTo asSafe( const TypeFrom& t )
{
  return ValueTypeConversionTraits<TypeTo,TypeFrom>::safeConvert(t);
}


template <class TypeTo> 
class asFunc {
  public:
  asFunc() {}

  template <class TypeFrom>
  inline TypeTo operator()(const TypeFrom &t) {
    return as<TypeTo>(t);
  }
};


//
// Standard specializations of ValueTypeConversionTraits
//


/** \brief Convert raw C string to std::string. */
template<int N>
class ValueTypeConversionTraits<std::string, char[N]> {
public:
  static std::string convert( const char t[] )
    { return std::string(t); }
  static std::string safeConvert( const char t[] )
    { return std::string(t); }
};

#ifdef HAVE_TEUCHOS_QD

/** \brief Convert qd_real to double. */
template <>
class ValueTypeConversionTraits<double, qd_real> {
public:
  inline static double convert( const qd_real t )
    { return to_double(t); }
  inline static double safeConvert( const qd_real t )
    { return to_double(t); }
};

/** \brief Convert qd_real to float. */
template <>
class ValueTypeConversionTraits<float, qd_real> {
public:
  inline static float convert( const qd_real t )
    { return (float)to_double(t); }
  inline static float safeConvert( const qd_real t )
    { return (float)to_double(t); }
};

/** \brief Convert qd_real to int. */
template <>
class ValueTypeConversionTraits<int, qd_real> {
public:
  inline static int convert( const qd_real t )
    { return to_int(t); }
  inline static int safeConvert( const qd_real t )
    { return to_int(t); }
};

/** \brief Convert qd_real to dd_real. */
template <>
class ValueTypeConversionTraits<dd_real, qd_real> {
public:
  inline static dd_real convert( const qd_real t )
    { return to_dd_real(t); }
  inline static dd_real safeConvert( const qd_real t )
    { return to_dd_real(t); }
};

/** \brief Convert dd_real to double. */
template <>
class ValueTypeConversionTraits<double, dd_real> {
public:
  inline static double convert( const dd_real t )
    { return to_double(t); }
  inline static double safeConvert( const dd_real t )
    { return to_double(t); }
};

/** \brief Convert dd_real to float. */
template <>
class ValueTypeConversionTraits<float, dd_real> {
public:
  inline static float convert( const dd_real t )
    { return (float)to_double(t); }
  inline static float safeConvert( const dd_real t )
    { return (float)to_double(t); }
};

/** \brief Convert dd_real to int. */
template <>
class ValueTypeConversionTraits<int, dd_real> {
public:
  inline static int convert( const dd_real t )
    { return to_int(t); }
  inline static int safeConvert( const dd_real t )
    { return to_int(t); }
};

#endif

// ToDo: Add more specializations as needed!


} // end namespace Teuchos


#endif // TEUCHOS_AS_HPP
