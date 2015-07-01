// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"

/*
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!     ATTENTION      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!! You might want to look at the BasicExample first. !!!!!!!!!!
 * !!!!!!!!!!!!! It's not ciritical, but it has some important     !!!!!!!!!!
 * !!!!!!!!!!!!! information that you might find helpful.          !!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

int main(int argc, char* argv[])
{
  /**
   * The end result of this example is exactly the same thing as we created
   * in the BasicExample. The difference is this time, instead of building
   * our GUI in source code, we're going to build it via XML. This approach
   * has a number of advantages.
   *    1. It's a lot cleaner. The source code approach can get pretty ugly
   *       especially when we get into some of the more complicated examples.
   *    2. Going along with the same vein of thought above, cleaner code is
   *       code that's easier to maintain.
   *    3. There's no need to recompile your GUI every time you change
   *       something. Since everything is loaded dynamically from the XML,
   *       there's no source code to recompile.
   */

  /**
   * The first thing we've gotta do is read in the ParameterList from the XML
   * so we can pass it to Optika. It's pretty easy to do.
   */
  Teuchos::RCP<Teuchos::ParameterList> myParameters;

  /**
   * Then we just call getInput! There's a little more to it, so let's
   * head on over to the inputs.xml file to see what's going on there.
   */
  Optika::getInput("inputs.xml", myParameters);

  return 0;
}

