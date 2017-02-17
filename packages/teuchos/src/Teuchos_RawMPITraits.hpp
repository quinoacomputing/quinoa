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

// ////////////////////////////////////////////////////////////////////////
// Teuchos_RawMPITraits.hpp

#ifndef TEUCHOS_RAW_MPI_TRAITS_H
#define TEUCHOS_RAW_MPI_TRAITS_H

#include "Teuchos_ConfigDefs.hpp"

/** \file Teuchos_RawMPITraits.hpp
 *  \brief Declaration of a templated traits class that returns raw MPI data types.
 */

namespace Teuchos {

/** \brief Templated traits class that allows a datatype to be used with MPI
 * that MPI can directly handle.
 *
 * A specialization of this traits class should only be created for datatypes
 * that can be directly handled by MPI in some way.  Note that this traits
 * class assumes that the datatype <tt>T</tt> is directly composed of
 * datatypes that MPI can directly handle.  This traits interface allows for
 * specializations to create user-defined <tt>MPI_Datatype</tt> and
 * <tt>MPI_Op</tt> objects to be returned from their static functions.
 *
 * \note 
 * <ul> <li> This class should not compile if it is instantiated by accident.  
 * 	<li> It should only be included if MPI is available and the MPI header 
 *		<b>must</b> be included before this header file.
 *	<li> Template specializations exist for datatypes: <tt>char</tt>, <tt>int</tt>,
 *		<tt>float</tt>, and <tt>double</tt>.
 *	<li> A partial template specialization exists for <tt>std::complex<T></tt>
 *    where it is assumed that the real type <tt>T</tt> is directly handlable
 *   with MPI.
 *  <li> Only sum reductions are supported for all data types.
 *  <li> The reductions max and min are only supported by datatypes whee
 *    <tt>ScalarTraits<T>::isComparable==true</tt> which is a compile-time
 *    boolean that can be used in template metaprogramming techniques.
 * </ul>
 */
template <class T> class RawMPITraits {
public:
  /** \brief Return the adjusted std::cout of items. */
  static int adjustCount(const int count) { bool *junk1; T *junk2 = &junk1; return 0; } // Should not compile!
	/** \brief Return the raw MPI data type of the template argument. */
	static MPI_Datatype type() { bool *junk1; T *junk2 = &junk1; return MPI_DATATYPE_NULL; } // Should not compile!
	/** \brief Return the MPI_Op object for a sum reduction */
	static MPI_Op sumOp() { bool *junk1; T *junk2 = &junk1; return MPI_OP_NULL; } // Should not compile!
	/** \brief Return the MPI_Op object for a max reduction */
	static MPI_Op maxOp() { bool *junk1; T *junk2 = &junk1; return MPI_OP_NULL; } // Should not compile!
	/** \brief Return the MPI_Op object for a min reduction */
	static MPI_Op minOp() { bool *junk1; T *junk2 = &junk1; return MPI_OP_NULL; } // Should not compile!
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Specialization of <tt>RawMPITraits</tt> for <tt>char</tt>
 */
template <> class RawMPITraits<char> {
public:
  /** \brief . */
  static int adjustCount(const int count) { return count; }
	/** \brief . */
	static MPI_Datatype type() { return MPI_CHAR; }
	/** \brief . */
	static MPI_Op sumOp() { return MPI_SUM; }
	/** \brief . */
	static MPI_Op maxOp() { return MPI_MAX; }
	/** \brief . */
	static MPI_Op minOp() { return MPI_MIN; }
};

/** \brief Specialization of <tt>RawMPITraits</tt> for <tt>int</tt>
 */
template <> class RawMPITraits<int> {
public:
  /** \brief . */
  static int adjustCount(const int count) { return count; }
	/** \brief . */
	static MPI_Datatype type() { return MPI_INT; }
	/** \brief . */
	static MPI_Op sumOp() { return MPI_SUM; }
	/** \brief . */
	static MPI_Op maxOp() { return MPI_MAX; }
	/** \brief . */
	static MPI_Op minOp() { return MPI_MIN; }
};

/** \brief Specialization of <tt>RawMPITraits</tt> for <tt>float</tt>
 */
template <> class RawMPITraits<float> {
public:
  /** \brief . */
  static int adjustCount(const int count) { return count; }
	/** \brief . */
	static MPI_Datatype type() { return MPI_FLOAT; }
	/** \brief . */
	static MPI_Op sumOp() { return MPI_SUM; }
	/** \brief . */
	static MPI_Op maxOp() { return MPI_MAX; }
	/** \brief . */
	static MPI_Op minOp() { return MPI_MIN; }
};

/** \brief Specialization of <tt>RawMPITraits</tt> for <tt>double</tt>
 */
template <> class RawMPITraits<double> {
public:
  /** \brief . */
  static int adjustCount(const int count) { return count; }
	/** \brief . */
	static MPI_Datatype type() { return MPI_DOUBLE; }
	/** \brief . */
	static MPI_Op sumOp() { return MPI_SUM; }
	/** \brief . */
	static MPI_Op maxOp() { return MPI_MAX; }
	/** \brief . */
	static MPI_Op minOp() { return MPI_MIN; }
};

/** \brief Partial specialization of <tt>RawMPITraits</tt> for <tt>std::complex<T></tt>.
 *
 * Note, <tt>maxOp()</tt> and <tt>minOp()</tt> are not supported by std::complex
 * numbers.
 *
 * ToDo: If a platform is found where this simple implementation does not work
 * then something else will have to be considered.
 */
template <class T> class RawMPITraits< std::complex<T> > {
public:
  /** \brief . */
  static int adjustCount(const int count) { return (2*count); }
	/** \brief . */
	static MPI_Datatype type() { return RawMPITraits<T>::type(); }
	/** \brief . */
	static MPI_Op sumOp() { return MPI_SUM; }
};

#endif // DOXYGEN_SHOULD_SKIP_THIS

} // namespace Teuchos

#endif // TEUCHOS_RAW_MPI_TRAITS_H
