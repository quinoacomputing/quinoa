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

#ifndef VECTOR_SPACE_TESTER_H
#define VECTOR_SPACE_TESTER_H

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Testing class for \c VectorSpace, \c Vector and \c VectorMutable.
 *
 * The purpose of this class is to test a \c VectorSpace object and the
 * \c VectorMutable objects that it creates.  The testing function
 * \c check_vector_space() calls all of the methods defined in the interfaces
 * \c %VectorSpace, \c %Vector and \c %VectorMutable and checks
 * many of the post conditions but not all.  It would be very difficult to 
 * completely verify every postcondition in every situation. 
 *
 * The behavior of the testing function check_vector_space() is strongly influenced
 * by a set of options (see \c VectorSpaceTester()).
 *
 * When writting new vector implementations, a developer is likely to spend a lot
 * of time debuggin while in this testing function.
 */
class VectorSpaceTester {
public:

  /// Members for option \c print_all_tests() (see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_all_tests );
#ifdef DOXYGEN_COMPILE
    ;
#endif		
  /// Members for option \c print_vectors() (see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_vectors );
#ifdef DOXYGEN_COMPILE
    ;
#endif		
  /// Members for option \c throw_exception() (see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception );
#ifdef DOXYGEN_COMPILE
    ;
#endif		
  /// Members for option \c num_random_tests() (see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_random_tests );
#ifdef DOXYGEN_COMPILE
    ;
#endif		
  /// Members for option \c () warning_tol(see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol );
#ifdef DOXYGEN_COMPILE
    ;
#endif		
  /// Members for option \c error_tol() (see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol );
#ifdef DOXYGEN_COMPILE
    ;
#endif		

  /** \brief Constructor (set default options).
   *
   * These default options are appropriate for even the largest vector spaces.
   */
  VectorSpaceTester(
    bool         print_all_tests  = false
    ,bool        print_vectors    = false
    ,bool        throw_exception  = true
    ,size_type   num_random_tests = 4
    ,value_type  warning_tol      = 1e-14
    ,value_type  error_tol        = 1e-10
    );

  /** \brief . */
  virtual ~VectorSpaceTester() {}

  /** \brief Run a vector space and the vectors it creates through a set of comprehensive tets.
   *
   * @param  space  [in] The vector space object to test.
   * @param  out    [in/out] If <tt>out != NULL</tt> then output will be sent to <tt>*out</tt>.
   *
   * The behavior of this function greatly depends on a number of options (see \c VectorSpaceTester()
   * for the default values for these options).  Access functions to set these options are provided
   * by the prototypes of the macro <tt>STANDARD_MEMBER_COMPOSITION_MEMBERS()</tt>.
   * <ul>
   * <li> <b><tt>print_all_tests(bool)</tt></b>:  If <tt>print_all_tests() == true</tt>, then some output will be sent to
   *      <tt>*out</tt> for every test performed.  This is useful to see all of tests that are performed and
   *      in debugging.
   * <li> <b><tt>print_vectors(bool)</tt></b>:  If <tt>print_vectors() == true</tt>, then all of the vectors will be printed
   *      that are created durring the tests.  This option is really only needed durring initial debugging
   *      and should only be used with small vector spaces since it will produce a lot of <tt>O(space.dim())</tt>
   *      output.
   * <li> <b><tt>throw_exception(bool)</tt></b>:  If <tt>throw_exception() == true</tt>, then any object that throws
   *      an unexpected exception will cause that exception to be thrown clear of of this function.  If
   *      <tt>out != NULL</tt> then the <tt>what()</tt> string will be printed to <tt>*out</tt> before the exception
   *      is rethrown.  If <tt>throw_exception() == false</tt>, then all exceptions will be caught, printed to 
   *      <tt>*out</tt> and then <tt>false</tt> is returned from the function.
   * <li> <b><tt>num_random_tests(int)</tt></b>:  This is the number of random tests to perform per category of test.
   *      A higher number will result is better validation but will consume more CPU time.
   * <li> <b><tt>warning_tol(value_type)</tt></b>:  Any test with a relative error greater than <tt>warning_tol()</tt> will
   *      result in a warning message printed to <tt>*out</tt>.
   * <li> <b><tt>error_tol(value_type)</tt></b>:  Any test with a relative error greater than <tt>erfor_tol()</tt> will
   *      result in an error message printed to <tt>*out</tt> and the function will immediatly return <tt>false</tt>.
   * </ul>
   */
  virtual bool check_vector_space(
    const VectorSpace &space
    ,std::ostream     *out
    ) const;

private:

  /** \brief . */
  void check_test(value_type err, std::ostream* out, bool* success) const;

}; // end class VectorSpaceTester

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_TESTER_H
