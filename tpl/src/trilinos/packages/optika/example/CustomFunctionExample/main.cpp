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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
/**
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 * !!!!!              ATTENTION              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!  PLEASE VIEW THE BASIC EXAMPLE FIRST BEFORE READING THIS    !!!!!!!!
 * !!!!  EXAMPLE. IT PROVIDES FUNDAMENTAL KNOWLEDGE THAT WILL BE    !!!!!!!!
 * !!!!  VERY HELPFUL IN UNDERSTANDING THIS EXAMPLE                 !!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */ 

/**
 * Don't worry about this little guy right now. We'll come back to him later.
 * We just needed to prototype it.
 */
void customFunction(Teuchos::RCP<const Teuchos::ParameterList> currentParams);

int main(int argc, char* argv[])
{
  /**
   * Sometimes, you would rather have a slightly different workflow than the 
   * one Optika offers by default. This is why Optika also offers a second,
   * tighter workflow that goes like this:
   * 	1. Construct a Parameter List of inputs.
   * 	2. Create a function which will run given the inputs.
   * 	3. Call getInput() like normal, except also pass the memory address of 
   * 	   the function.
   * 	4. Everytime the user clicks the action button on the GUI, your custom 
   * 	   function will run.
   * 	5. The user may click quit when they are finished and control will be 
   * 	   returned to you.
   * This alternative workflow allows the user to quickly tweak and re-try 
   * different parameter configurations.
   */

  /**
   * We start off like normal. I'm actually just going to repeat what we did 
   * in the BasicExample
   */
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::StringValidator;
  using Teuchos::tuple;
  using Teuchos::Array;
  using Teuchos::rcp;
  RCP<ParameterList> My_List = rcp(new ParameterList);

  My_List->set(
    "Max Iters", 
    1550, 
    "Determines the maximum number of iterations in the solver");
  My_List->set(
    "Tolerance", 
    1e-10, 
    "The tolerance used for the convergence check");
  
  RCP<StringValidator> solverValidator = 
     rcp(new StringValidator(tuple<std::string>("GMRES", "CG", "TFQMR")));
  My_List->set( 
    "Solver", 
    "GMRES", 
    "The type of solver to use.", 
    solverValidator);

  Array<double> doubleArray( 10, 0.0 );
  My_List->set(
    "Initial Guess", 
    doubleArray, 
    "The initial guess as a RCP to an array object.");

  ParameterList& Prec_List = My_List->sublist(
    "Preconditioner",
    false,
    "Sublist that defines the preconditioner.");

  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set(
    "Drop Tolerance",
    1e-3,
    "The tolerance below which entries from the "
    "factorization are left out of the factors.");

  /**
   * Here is where things get switched up a bit. Along with the list, we pass 
   * along the address of the "customFunction" (the Teuchos::null is there because we
   * don't want to specify any dependencies for this GUI). Let's skip down a bit and take 
   * a look at the function.
   */
  Optika::getInput(My_List, Teuchos::null, &customFunction);

  RCP<Teuchos::FancyOStream> out = 
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::writeParameterListToXmlOStream(*My_List, *out);
  
  return 0;
}

/**
 * The custom function must always have this signature (returning void and 
 * taking a single RCP<const ParameterList> parameter). When the user clicks 
 * action button the function will be called. The current values for the 
 * parameter list are what will be given as the argument. Here we simply read 
 * what some of the settings the user has selected and print them out. But in 
 * theory, you could do anything you want with the parameters at this point.
 * The GUI will terminate once the user exits (by some means of closing the GUI
 * window), and control will be returned back to what ever function called 
 * getInput() (in this case the main()).
 */
void customFunction(Teuchos::RCP<const Teuchos::ParameterList> currentParams){
  std::cout << "The Solver Choosen was: " << 
    currentParams->get<std::string>("Solver") << std::endl;

  std::cout << "The Tolerence choosen was: " << 
    currentParams->get<double>("Tolerance") << std::endl; 

  std::cout << "The Preconditioner type was set to: " 
    << currentParams->sublist("Preconditioner").get<std::string>("Type") << 
    std::endl;
}

