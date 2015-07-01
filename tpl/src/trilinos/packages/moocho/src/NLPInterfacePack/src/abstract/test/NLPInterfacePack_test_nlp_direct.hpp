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

#ifndef TEST_NLP_FIRST_ORDER_DIRECT_H
#define TEST_NLP_FIRST_ORDER_DIRECT_H

#include <iosfwd>

#include "NLPInterfacePack_Types.hpp"

namespace OptionsFromStreamPack {
  class OptionsFromStream;
}

namespace NLPInterfacePack {

/** \brief Test an <tt>NLPDirect</tt> object.
 *
 * @param  nlp  [in/out] %NLP object being tested.
 * @param  options
 *              [in] If <tt>options != NULL</tt> then the options to use are extracted
 *              from <tt>*options</tt>.  If <tt>options == NULL</tt> then a default set
 *              of options will be used that will be appropriate for even the largest %NLP
 *              (see below).
 * @param  out  [in/out] If <tt>out != NULL</tt> then output will be set to <tt>*out</tt>.
 *              The amount of output sent to <tt>*out</tt> depends on the options selected.
 *              If <tt>out == NULL</tt> then no output is produced.
 *
 * This function uses the testing classes <tt>\ref AbstractLinAlgPack::VectorSpaceTester "VectorSpaceTester"</tt>
 * <tt>\ref NLPInterfacePack::NLPTester "NLPTester"</tt> and
 * <tt>\ref NLPInterfacePack::NLPDirectTester "NLPDirectTester"</tt> to perform many thorough tests
 * of an input <tt>\ref NLPInterfacePack::NLPDirect "NLPDirect"</tt> object.
 * The vector spaces exposed by <tt>\ref NLPInterfacePack::NLP "NLP"</tt> are thoroughly tested by the <tt>VectorSpaceTester</tt>
 * class.
 *
 * The options groups "VectorSpaceTester" (see <tt>\ref AbstractLinAlgPack::VectorSpaceTesterSetOptions "VectorSpaceTesterSetOptions"</tt>),
 * "%NLPTester" (see <tt>\ref NLPInterfacePack::NLPTesterSetOptions "NLPTesterSetOptions"</tt>), "%CalcFiniteDiffProd"
 * (see <tt>\ref NLPInterfacePack::CalcFiniteDiffProdSetOptions "CalcFiniteDiffProdSetOptions"</tt>) and "%NLPDirectTester"
 * (see <tt>\ref NLPInterfacePack::NLPDirectTesterSetOptions "NLPDirectTesterSetOptions"</tt>) are looked for in
 * in <tt>*options</tt> (if <tt>options != NULL</tt>) in order to extract options to use for this testing function and the other testing
 * objects.
 */
bool test_nlp_direct(
  NLPDirect                                     *nlp
  ,OptionsFromStreamPack::OptionsFromStream     *options
  ,std::ostream                                 *out
  );

} // end namespace NLPInterfacePack

#endif // TEST_NLP_FIRST_ORDER_DIRECT_H
