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

#ifndef NLP_INTERFACE_PACK_NLP_TESTER_H
#define NLP_INTERFACE_PACK_NLP_TESTER_H

#include <iosfwd>

#include "NLPInterfacePack_Types.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace NLPInterfacePack {

/** \brief Testing class for base NLP interface.
 *
 * This class is little more than a unit tester for the
 * <tt>NLP</tt> base interface.  This class will call all
 * of the <tt>%NLP</tt> methods and print out quanities if asked to.
 * This class simply validates the pre and post conditions
 * for all of the methods.  In that this class is useful.
 */
class NLPTester {
public:

  /// Members for option \c print_all() (see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_all );
#ifdef DOXYGEN_COMPILE
    ;
#endif		
  /// Members for option \c throw_exception() (see Teuchos_StandardMemberCompositionMacros.hpp).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception );
#ifdef DOXYGEN_COMPILE
    ;
#endif		
  ///	Constructor (default options)
  NLPTester(
    bool     print_all        = false
    ,bool    throw_exception  = true
    );
 
  /** \brief Test the NLP interface as the given base point xo.
   *
   * @param  nlp  [in/out] The <tt>NLP</tt> object being tested.  The pointers returned by
   *              <tt>nlp->get_f()</tt>, <tt>nlp->get_c()</tt> and <tt>nlp->get_h()</tt> will
   *              be preserved on output and will not be modified by this function..  The
   *              %NLP must be initialized before input.
   * @param  xo   [in] Base point for the unknown variables to test the calcuation methods at.
   * @param  print_all_warnings
   *              [in] Determines if warnings for all of the comparison tests are printed or not.
   *              Warning: may cause as much as <i>O(</i><tt>this->n()<tt><i>)</i> output.
   * @param  out  [in/out] If <tt>out != NULL</tt> any and all output will be sent here.  If
   *              <tt>out == NULL</tt> then no output will be produced.
   *
   * @return Returns \c true if all of the tests checked out and no unexpected exceptions were
   * thrown.
   *
   * The behavior of this method depends on the options \c print_all() and
   * \c throw_exception() and the input arguments \c print_all_warnings and
   * \c out.
   * <ul>
   * <li> <b><tt>throw_exception(bool)</tt></b>:
   * If <tt>throw_exception()</tt> == true</tt>, then if any of the objects within
   * this function throw exceptions, these exceptions will be be thrown clean
   * out of this function for the caller to handle.  If <tt>throw_exception()</tt> == false</tt>,
   * then if any object throws an exception, the exception is caught and this this function will
   * return <tt>false</tt>.  In any case an error message will be printed
   * to <tt>*out</tt> (if <tt>out != NULL</tt) before leaving the function (by \c return or \c throw).
   * <li> <b><tt>print_all(bool)</tt></b>:
   * If <tt>print_all() == true</tt> then all of the computed quantities will but dumped to \c out.
   * Note that this is a useful option for initial debugging of small NLPs but not a good idea for
   * larger NLPs as it will result in an excessive amount of output.
   * </ul>
   */
  bool test_interface(
    NLP                     *nlp
    ,const Vector           &xo
    ,bool                   print_all_warnings
    ,std::ostream           *out
    );

}; // end class NLPTester

} // end namespace NLPInterfacePack

#endif // NLP_INTERFACE_PACK_NLP_TESTER_H
