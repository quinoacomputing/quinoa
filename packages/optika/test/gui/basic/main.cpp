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
#include "Teuchos_FancyOStream.hpp"

namespace Optika{


void RUN_BASIC_OPTIKA_TESTS(){

  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::StringToIntegralParameterEntryValidator;
  using Teuchos::tuple;
  using Teuchos::rcp;

  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

  //Basic Test
  RCP<ParameterList> My_List = RCP<ParameterList>(new ParameterList);

  double *pointer = 0;
  My_List->set("Double pointer", pointer);
  My_List->set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
  My_List->set("Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  RCP<StringToIntegralParameterEntryValidator<int> >
    solverValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "GMRES", "CG", "TFQMR" )
        ,"Solver"
        )
      );
  My_List->set(
    "Solver"
    ,"GMRES" // This will be validated by solverValidator right here!
    ,"The type of solver to use."
    ,solverValidator
    );

   RCP<EnhancedNumberValidator<int> > awesomenessValidator = 
   RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>(0,10));
    My_List->set("Awesomeness", 5, "Rate the awesomeness!!!", awesomenessValidator);

  Array<double> testArray( 10, 0.0 );
  
  My_List->set("Initial Guess", testArray, "The initial guess as a RCP to an array object.");

  ParameterList&
    Prec_List = My_List->sublist("Preconditioner",false,"Sublist that defines the preconditioner.");

  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.");

  getInput(My_List);

  writeParameterListToXmlOStream(*My_List, *out);


}


}

int main(){
  Optika::RUN_BASIC_OPTIKA_TESTS();
  return 0;
}

