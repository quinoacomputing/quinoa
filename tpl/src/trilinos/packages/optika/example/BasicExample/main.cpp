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
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::tuple;
using Teuchos::StringValidator;
using Teuchos::Array;
using Teuchos::VerboseObjectBase;
using Teuchos::FancyOStream;

int main(int argc, char* argv[])
{
 /*
  * Welcome to the Optika Package!
  *
  * This package was designed to assist in the rapid development of GUIs for 
  * existing and new projects using the Trilinos Framework. Using the 
  * ParameterList class found in the Teuchos package and the DependencySheet
  * class also provided by the Teuchos Package, Optika will allow you to use 
  * ParameterLists to define a set of values you wish to obtain from
  * the user. You may then pass this ParameterList to the function getInput. 
  * This function will dynamically generate a GUI based on your ParameterList, 
  * display the GUI to the user, obtain input from the user, and then store 
  * the users input back into the ParameterList. 
  *
  * Let's take a look at an example to see how this all works.
  *
  * Before you Start:
  * We recommend you have at least a basic understanding of what a RCP is. 
  * While not * crucial to the understanding of these examples, undestanding 
  * RCPs allow you to more easily read what is going on in the examples.
  */


  /* 
   * First we create an empty parameter list. We will use this to define
   * all of the parameters we wish to obtain from the user. This type of 
   * ParameterList is commonly known as the "Valid Parameter List".
   */
  RCP<ParameterList> My_List = RCP<ParameterList>(new ParameterList);

  /* 
   * Creating parameters in this list can be easily done using the set function.
   * The first argument is the name of the parameter. The second is the default
   * value for the parameter. The third is a short description of what the 
   * parameter is for.
   */
  My_List->set(
    "Max Iters", 
    1550, 
    "Determines the maximum number of iterations in the solver");

  My_List->set(
    "Tolerance", 
    1e-10,
    "The tolerance used for the convergence check");
  
  /* 
   * Validators are useful for restricting the set of values that may be used 
   * for a given parameter. For the "Solver" option, we will create a 
   * validator. Here we use a StringValidator and a tuple to specify which 
   * string values are valid for the "Solver" option.
   */
  RCP<StringValidator> solverValidator = 
     RCP<StringValidator>(new StringValidator(
      tuple<std::string>("GMRES", "CG", "TFQMR")));

  My_List->set(
    "Solver", 
    "GMRES", 
    "The type of solver to use.", 
    solverValidator);

  /* 
   * The Optika Package can also handle Teuchos Arrays.
   * Here we create a Array object of 10 doubles
   * representing an initial guess for a linear solver.
   */
  Array<double> doubleArray( 10, 0.0 );
  My_List->set(
    "Initial Guess", 
    doubleArray, 
    "The initial guess as an array object.");

  /* 
   * We can also create a hieiarchy of parameters by using sublists. Here we 
   * create a sublist called Prec_List. Prec_List will be contained within 
   * My_List.
   */
  ParameterList& Prec_List = My_List->sublist(
    "Preconditioner",
    false,
    "Sublist that defines the preconditioner.");

  /*
   * Now this Prec_List can be filled with other parameters:
   */
  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n"
                "factorization are left out of the factors.");

  /*
   * The getInput function starts up an Optika GUI and lets the user start to 
   * input parameter values. When the user  has completed their data entry, 
   * the function will finish right after all of the input values are stored 
   * in My_List.
   */
  Optika::getInput(My_List);

  /*
   * Here we can print out what the user entered in nice XML format.
   */
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  writeParameterListToXmlOStream(*My_List, *out);
  

  /*
   * A Few Final Notes
   *
   * -After calling the getInput function, any parameter in My_List has the 
   *  potential to have been modified.
   *  That said, no new parameters or ParameterLists will have been added and 
   *  none will have been removed.
   *
   * -The GUI can only handle certain types of parameters. They are:
   *	int
   * 	short
   * 	double
   * 	float
   * 	bool
   * 	std::string
   * 	Array<int>
   * 	Array<short>
   * 	Array<double>
   * 	Array<float>
   * 	Array<string>
   * If you give it a ParameterList containing a parameter that is not of one 
   * of the types specified above, the parameter will still be displayed in 
   * the GUI. However, the user will not be able to modify it's value.
   *
   * That's it for now. Be sure check out the other examples to see some of 
   * the more advanced features of the Optika package. If you have any 
   * suggestions or feature requests, please send them to klnusbaum@gmail.com.
   */
  Teuchos::ParameterEntryXMLConverterDB::printKnownConverters(*out);
  return 0;
}

