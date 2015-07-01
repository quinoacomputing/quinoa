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

#ifndef DECOMPOSITION_SYSTEM_TESTER_H
#define DECOMPOSITION_SYSTEM_TESTER_H

#include <iosfwd>

#include "ConstrainedOptPack_Types.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Testing class for \c DecompositionSystem interface.
 *
 * This testing class is basically a unit tester for \c DecompositionSystem.  The method \c test_decomp_system()
 * runs many different tests to validate the interface and the matrix objects associated with the interface.
 * The method \c test_decomp_system() should only be called after
 * <tt>decomp_sys\ref DecompositionSystem::update_decomp ".update_decomp(...)"</tt> is called on the
 * <tt>DecompositionSystem</tt> object <tt>decomp_sys</tt>.  The output decomposition matrices are passed
 * through a series of tests.  The compatibility of the matrices described in the postconditions for
 * <tt>DecompositionSytem::update_decomp()</tt> are also checked in a series of tests.  If the method
 * \c test_decomp_system() returns \c true, then the client can feel fairly confident that the
 * things are functioning properly (but this is not guaranteed of course).
 *
 * The tests performed by this testing class are designed to allow some validation for even the largest systems
 * and will produce various levels of output so as to be usefull in debugging.
 *
 * ToDo:  Finish documentation!
 */
class DecompositionSystemTester {
public:

  /** @name Set and access options */
  //@{

  /** \brief . */
  enum EPrintTestLevel {
    PRINT_NOT_SELECTED =0  ///< The print option has not been selected (will default to PRINT_NONE if not set)
    ,PRINT_NONE        =1  ///< Don't print anything
    ,PRINT_BASIC       =2  ///< Print only very basic info
    ,PRINT_MORE        =3  ///< Print greater detail about the tests.
    ,PRINT_ALL         =4  ///< Print everything all the tests in great detail but output is independent of problem size.
  };
  /// Set the level of output produced durring tests.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EPrintTestLevel, print_tests );
  /// Set whether matrices, vectors ect. are printed (warning, this may be a lot of output for larger systems).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );
  /// Set whether an exception that is thrown is thrown clear out of the testing function or not.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception );
  /// Set the number of random test cases created.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_random_tests );
  /// Set the relative tolerance for numerical tests of matrix-vector multiplication above which to print a warning.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, mult_warning_tol );
  /// Set the relative tolerance for numerical tests  of matrix-vector multiplication above which to return \c false from the testing function.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, mult_error_tol );
  /// Set the relative tolerance for numerical tests of matrix-vector solves above which to print a warning.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, solve_warning_tol );
  /// Set the relative tolerance for numerical tests  of matrix-vector solves above which to return \c false from the testing function.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, solve_error_tol );

  //@}

  ///	Constructor (default options)
  DecompositionSystemTester(
    EPrintTestLevel  print_tests       = PRINT_NOT_SELECTED
    ,bool            dump_all          = false
    ,bool            throw_exception   = true
    ,size_type       num_random_tests  = 1
    ,value_type      mult_warning_tol  = 1e-14
    ,value_type      mult_error_tol    = 1e-8
    ,value_type      solve_warning_tol = 1e-14
    ,value_type      solve_error_tol   = 1e-8
    );
 
  /** \brief Test a \c DecompositionSystem object after <tt>DecompositionSystem::update_basis()</tt> is called.
   *
   * @param  decomp_sys
   *              [in] The \c DecompositionSystem object that \c DecompositionSystem::update_basis() was called on.
   * @param  Gc   [in] Must be matrix that was passed to <tt>decomp_sys.update_decomp(Gc,...)</tt>.
   * @param  Z    [in] Must be matrix that was passed to and returned from
   *              <tt>decomp_sys.update_decomp(...,Z,...)</tt>.
   * @param  Y    [in] Must be matrix that was passed to and returned from
   *              <tt>decomp_sys.update_decomp(...,Y,...)</tt>.
   * @param  R    [in] Must be matrix that was passed to and returned from
   *              <tt>decomp_sys.update_decomp(...,R,...)</tt>.
   * @param  Uz   [in] If <tt>Uz!=NULL</tt> then Must be matrix that was passed to and returned from
   *              <tt>decomp_sys.update_decomp(...,Uz,...)</tt>.
   * @param  Uy   [in] If <tt>Uy!=NULL</tt> then Must be matrix that was passed to and returned from
   *              <tt>decomp_sys.update_decomp(...,Uy,...)</tt>.
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
   *      out of this function for the caller to handle.
   *      If <tt>throw_exception()</tt> == false</tt>, then if any object throws an exception,
   *      the exception is caught and this this function will return <tt>false</tt>.
   *      In any case an error message will be printed to <tt>*out</tt> (if <tt>out != NULL</tt)
   *      before leaving the function (by \c return or \c throw).
   * <li> <b><tt>dump_all(bool)</tt></b>:
   *      If <tt>dump_all() == true</tt> then all of the computed quantities will but dumped to \c out.
   *      Note that this is a useful option for initial debugging of small systems but not a good idea for
   *      larger systems as it will result in an excessive amount of output.
   * <li> ToDo: Add rest of options!
   * </ul>
   */
  bool test_decomp_system(
    const DecompositionSystem   &decomp_sys
    ,const MatrixOp             &Gc
    ,const MatrixOp             *Z
    ,const MatrixOp             *Y
    ,const MatrixOpNonsing      *R
    ,const MatrixOp             *Uz
    ,const MatrixOp             *Uy
    ,std::ostream               *out
    );

}; // end class DecompositionSystemTester

} // end namespace ConstrainedOptPack

#endif // DECOMPOSITION_SYSTEM_TESTER_H
