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

#ifndef BASIS_SYSTEM_TESTER_H
#define BASIS_SYSTEM_TESTER_H

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Testing class for \c BasisSystem interface.
 *
 * This testing class is basically a unit tester for \c BasisSystem.  The method \c test_basis_system()
 * runs many different tests to validate the interface and the objects allocated with the interface.
 * The method \c test_basis_system() should only be called after
 * <tt>basis_sys\ref BasisSystem::update_basis ".update_basis(...)"</tt> is called on the <tt>BasisSystem</tt>
 * object <tt>basis_sys</tt>.  The output basis matrix \a C and/or direct sensitivity matrix \a D are passed through
 * a series of tests using the testing classes <tt>MatrixOpNonsingTester</tt> and <tt>MatrixWithOpTester</tt>
 * respectively.  The compatibility of the matrices \c Gc, \c C and/or \c D are also checked in a series of
 * tests.  If the method \c test_basis_system() returns true, then the client can feel fairly confident that the
 * basis matrix object is functioning properly.
 *
 * The tests performed by this testing class are designed to allow some validation for even the larges systems
 * and will produce various levels of output so as to be usefull in debugging.
 *
 * ToDo:  Finish documentation!
 */
class BasisSystemTester {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum EPrintTestLevel {
    PRINT_NOT_SELECTED = 0  ///< The print option has not been selected (will default to PRINT_NONE if not set)
    ,PRINT_NONE        = 1  ///< Don't print anything
    ,PRINT_BASIC       = 2  ///< Print only very basic info
    ,PRINT_MORE        = 3  ///< Print greater detail about the tests.
    ,PRINT_ALL         = 4  ///< Print everything all the tests in great detail but output is independent of problem size.
  };

  //@}

  /** @name Set and access options */
  //@{

  /// Set the level of output produced durring tests.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EPrintTestLevel, print_tests );
  /// Set whether matrices, vectors ect. are printed (warning, this may be a lot of output for larger systems).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );
  /// Set whether an exception that is thrown is thrown clear out of the testing function or not.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception );
  /// Set the number of random test cases created.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_random_tests );
  /// Set the relative tolerance for numerical tests above which to print a warning.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol );
  /// Set the relative tolerance for numerical tests above which to return false from the testing function.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol );

  //@}

  /** @name Constructors / initializers */
  //@{

  ///	Constructor (default options)
  BasisSystemTester(
    EPrintTestLevel  print_tests      = PRINT_NOT_SELECTED
    ,bool            dump_all         = false
    ,bool            throw_exception  = true
    ,size_type       num_random_tests = 1
    ,value_type      warning_tol      = 1e-14
    ,value_type      error_tol        = 1e-8
    );

  //@}

  /** @name Test basis system */
  //@{
 
  /** \brief Test a \c BasisSystem object after <tt>BasisSystem::update_basis()</tt> is called.
   *
   * @param  basis_sys
   *              [in] The \c BasisSystem object that \c BasisSystem::update_basis() was called on.
   * @param  Gc   [in] Matrix \c Gc that was passed into \c basis_sys.update_basis() (if not \c NULL).
   * @param  C    [in] Matrix \c C that was passed in and out of \c basis_sys.update_basis() (if not \c NULL).
   * @param  N    [in] If not \c NULL, then this must the matrix \a N described in the documentation for
   *              \c BasisSystem.  This allows a matrix object created by the client to be check out here
   *              also.
   * @param  D    [in] Matrix \c D that was passed in and out of \c basis_sys.update_basis() (if not \c NULL).
   *              Actually, this can be any matrix object that the client may want to define that takes the
   *              role of \c D.  Such a matrix object can be tested here in this function along with the rest
   *              of the matrices.
   * @param  GcUP [in] Matrix \c GcUP that was passed in and out of \c basis_sys.update_basis() (if not \c NULL).
   *              Actually, this can be any matrix object that the client may want to define that takes the
   *              role of \c GcUP.  Such a matrix object can be tested here in this function along with the rest
   *              of the matrices.
   * @param  out  [in/out] If <tt>out != NULL</tt> any and all output will be sent here.  If
   *              <tt>out == NULL</tt> then no output will be produced.
   *
   * @return Returns \c true if all of the tests checked out and no unexpected exceptions were
   * thrown.
   *
   * The behavior of this method depends on a set of options and the input arguments.
   * <ul>
   * <li> <b><tt>throw_exception(bool)</tt></b>:
   *      If <tt>throw_exception()</tt> == true</tt>, then if any of the objects within
   *      this function throw exceptions, these exceptions will be be thrown clean
   *      out of this function for the caller to handle.  If <tt>throw_exception()</tt> == false</tt>,
   *      then if any object throws an exception, the exception is caught and this this function will
   *      return <tt>false</tt>.  In any case an error message will be printed
   *      to <tt>*out</tt> (if <tt>out != NULL</tt) before leaving the function (by \c return or \c throw).
   * <li> <b><tt>dump_all(bool)</tt></b>:
   *      If <tt>dump_all() == true</tt> then all of the computed quantities will but dumped to \c out.
   *      Note that this is a useful option for initial debugging of small systems but not a good idea for
   *      larger systems as it will result in an excessive amount of output.
   * <li> ToDo: Add rest of options!
   * </ul>
   */
  bool test_basis_system(
    const BasisSystem           &basis_sys
    ,const MatrixOp             *Gc
    ,const MatrixOpNonsing      *C
    ,const MatrixOp             *N
    ,const MatrixOp             *D
    ,const MatrixOp             *GcUP
    ,std::ostream               *out
    );

  //@}

}; // end class BasisSystemTester

} // end namespace AbstractLinAlgPack

#endif // BASIS_SYSTEM_TESTER_H
