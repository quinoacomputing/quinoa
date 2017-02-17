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

#ifndef TEUCHOS_ARRAY_ARG_HPP
#define TEUCHOS_ARRAY_ARG_HPP


#include "Teuchos_Assert.hpp"


namespace Teuchos {


/** \defgroup Teuchos_Array_Arguments Utility functions for passing arrays into argument lists [Deprecated].

\brief The purpose of this utility is to make passing arrays into argument lists easier [Deprecated].

\ingroup Teuchos_Depreicated_grp

Declaring arrays outside of a function just to pass a (small) list of
values into a function can be tiresome. The templated function
<tt>arrayArg()</tt> simplifies this process.  With this function you
can construct (using stack memory not dynamically allocated memory) an
array of data to be passed into a function.

For example, consider the following function prototype:

\code
void f( const int x_size, const double x[] );
\endcode

which takes an array of <tt>double</tt>s of length <tt>x_size</tt>.
Generally, to call this function one would have to first declare an
array and then call the function as:

\code
void f()
{
  ...
  const double x[] = { 1.0, 2.0, 3.0 };
  f( 3, x );
  ...
\endcode

Now, however, one can create the array in the call to <tt>f()</tt> as:

\code
void f()
{
  ...
  f( 3, arrayArg<double>(1.0,2.0,3.0)() );
  ...
}
\endcode

In the above situation, one may be able to write the call as:

\code
void f()
{
  ...
  f( 3, arrayArg(1.0,2.0,3.0) );
  ...
}
\endcode

but the former, slightly more verbose, version is to be preferred
since it makes explicit what type of array is being created and insures
that the compiler will not get confused about the final (implicit) conversion
to a raw <tt>const double*</tt> pointer.

Note that a copy is made of the array arguments before they are passed into
the function so care must be taken when using <tt>arrayArg()</tt> to pass
a non-<tt>const</tt> input-output or output-only array of objects.  For example,
consider the following function:

\code
void f2( const int y_size, double y[] );
\endcode

The above function <tt>f2()</tt> modifies the objects in the array <tt>y[]</tt>.
If this function is attempted to be called as:

\code
void g2()
{
  double a, b, c;
  f2( 3, arrayArg(a,b,c)() );
}
\endcode

then the objects <tt>a</tt>, <tt>b</tt> and <tt>c</tt> will not be
modified as might be expected.  Instead, this function must be called as:

\code
void g2()
{
  double y[3];
  f2( 3, y );
  double a=y[0], b=y[1], c=y[2];
}
\endcode

However, the <tt>arrayArg()</tt> function can be used to pass
an array of pointers to non-<tt>const</tt> objects.  For example,
consider the function:

\code
void f3( const int y_size, double* y[] );
\endcode

which modifies an array of <tt>double</tt> objects through pointers.
We could then call this function as:

\code
void g3()
{
  double a, b, c;
  f2( 3, arrayArg(&a,&b,&c)() );
}
\endcode

which will result in objects <tt>a</tt>, <tt>b</tt> and <tt>c</tt>
being modified correctly.

Warning! Never try to pass an array of references (which should almost
never be used anyway) using <tt>arrayArg()</tt>.  This will result in the
copy constructor being called which is almost never a desirable situation.

The <tt>arrayArg()</tt> function is overloaded to accept 1, 2, 3, 4, 5 and 6
arguments.  If more elements are needed, then more overrides are easy to add.

*/

/** \brief Utility class that allows arrays to be passed into argument list.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<int N, class T>
class ArrayArg {
public:
  /// Basic constructor taking a copy of the \c array of length \c N
  ArrayArg( T array[] ) { std::copy( array, array+N, array_ ); }

  /// Return a \c const pointer to the internal array
  T* operator()() { return array_; }

  /// Return a \c const pointer to the internal array
  operator T* () { return array_; }

private:
  T array_[N]; //  Can't be a const array!
}; // class Array1DArg

/** \brief Return an array with 1 member.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline TEUCHOS_DEPRECATED ArrayArg<1,T> arrayArg( T t1 )
{
  T array[] = { t1 };
  return ArrayArg<1,T>(array);
}

/** \brief Return an array with 2 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline TEUCHOS_DEPRECATED ArrayArg<2,T> arrayArg( T t1, T t2 )
{
  T array[] = { t1, t2 };
  return ArrayArg<2,T>(array);
}

/** \brief Return an array with 3 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline TEUCHOS_DEPRECATED ArrayArg<3,T> arrayArg( T t1, T t2, T t3 )
{
  T array[] = { t1, t2, t3 };
  return ArrayArg<3,T>(array);
}

/** \brief Return an array with 4 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline TEUCHOS_DEPRECATED ArrayArg<4,T> arrayArg( T t1, T t2, T t3, T t4 )
{
  T array[] = { t1, t2, t3, t4 };
  return ArrayArg<4,T>(array);
}

/** \brief Return an array with 5 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline TEUCHOS_DEPRECATED ArrayArg<5,T> arrayArg( T t1, T t2, T t3, T t4, T t5 )
{
  T array[] = { t1, t2, t3, t4, t5 };
  return ArrayArg<5,T>(array);
}

/** \brief Return an array with 6 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline TEUCHOS_DEPRECATED ArrayArg<6,T> arrayArg( T t1, T t2, T t3, T t4, T t5, T t6 )
{
  T array[] = { t1, t2, t3, t4, t5, t6 };
  return ArrayArg<6,T>(array);
}

} // namespace Teuchos

#endif // TEUCHOS_ARRAY_ARG_HPP
