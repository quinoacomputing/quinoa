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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Optika_GUI.hpp"

void print(Teuchos::RCP<const Teuchos::ParameterList> theList){
  
  Teuchos::RCP<Teuchos::FancyOStream> out = 
    Teuchos::VerboseObjectBase::getDefaultOStream();
  writeParameterListToXmlOStream(*theList, *out);
}

namespace Optika{


void RUN_OPTIKA_DEPENDENCY_AND_EXEC_TEST(){
 using Teuchos::FancyOStream;
 using Teuchos::VerboseObjectBase;
 using Teuchos::StringToIntegralParameterEntryValidator;
 using Teuchos::StringValidatorDependency;
 using Teuchos::ArrayNumberValidator;
 using Teuchos::StringVisualDependency;
 using Teuchos::NumberVisualDependency;
 using Teuchos::NumberArrayLengthDependency;
 using Teuchos::BoolValidatorDependency;
 using Teuchos::BoolVisualDependency;
 using Teuchos::SubtractionFunction;
 using Teuchos::RangeValidatorDependency;
 using Teuchos::tuple;
 using Teuchos::rcp;


 RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
 RCP<ParameterList> My_deplist = rcp(new ParameterList);
 RCP<DependencySheet> depSheet1 = rcp(new DependencySheet);

  RCP<StringToIntegralParameterEntryValidator<int> >
    stringFoodTypeValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "Cheese", "Soda", "Chips" )
        ,"Food Type"
        )
      );

  RCP<StringToIntegralParameterEntryValidator<int> >
    cheeseValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "Swiss", "American", "Super Awesome Cheese" )
        ,"Food Selector"
        )
      );
  RCP<StringToIntegralParameterEntryValidator<int> >
    sodaValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "Pepsi", "Coke", "Kurtis Cola", "Bad Cola" )
        ,"Food Selector"
        )
      );
  RCP<StringToIntegralParameterEntryValidator<int> >
    chipsValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "Lays", "Doritos", "Kurtis Super Awesome Brand" )
        ,"Food Selector"
        )
      );

  StringValidatorDependency::ValueToValidatorMap testValidatorMap1;
  testValidatorMap1["Cheese"] = cheeseValidator;
  testValidatorMap1["Soda"] = sodaValidator;
  testValidatorMap1["Chips"] = chipsValidator;


  ParameterList&
    stringValiDepList = My_deplist->sublist("String Validator Dependency", false, "String Validator Dependency testing list.\nWorking June 27th 2009");
  stringValiDepList.set("Food Selector", "Swiss", "select the food you want", cheeseValidator);
  stringValiDepList.set("Food Type", "Cheese", "String Validator Dependency Tester", stringFoodTypeValidator);
  depSheet1->addDependency(RCP<StringValidatorDependency>(
  	new StringValidatorDependency(
    My_deplist->getEntryRCP("Food Type"),
	  My_deplist->getEntryRCP("Food Selector"),
	  testValidatorMap1 )));
 

  RCP<StringToIntegralParameterEntryValidator<int> >
    stringRangeValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "1-10", "10-33", "50-60" )
        ,"Range selector"
        )
      );
  ParameterList&
    stringValiDepList2 = My_deplist->sublist("String Validator Dependency (other validators)", false, "String Validator Dependency testing list for EnhancedNumber Validators.");
  stringValiDepList2.set("Range selector", "1-10", "selects the range to validate", stringRangeValidator);
  RCP<EnhancedNumberValidator<int> > range110Vali = 
  	rcp(new EnhancedNumberValidator<int>(1,10));
  RCP<EnhancedNumberValidator<int> > range1033Vali = 
  	rcp(new EnhancedNumberValidator<int>(10,33));
  RCP<EnhancedNumberValidator<int> > range5060Vali = 
  	rcp(new EnhancedNumberValidator<int>(50,60));
  StringValidatorDependency::ValueToValidatorMap rangeValidatorMap1;
  rangeValidatorMap1["1-10"] = range110Vali;
  rangeValidatorMap1["10-33"] = range1033Vali;
  rangeValidatorMap1["50-60"] = range5060Vali;
  stringValiDepList2.set("RangeValue", 3, "the value of the range", range110Vali);
  depSheet1->addDependency(RCP<StringValidatorDependency>(
  	new StringValidatorDependency(
	My_deplist->getEntryRCP("Range selector"),
	My_deplist->getEntryRCP("RangeValue"),
	rangeValidatorMap1)));


  ParameterList&
    boolValidatorDepList = My_deplist->sublist("Bool Validator Dependency List", false, "Bool Validator Dependency testing list.\nWorking June 27th 2009");
  boolValidatorDepList.set("Use Validator?", true, "truns the validator on and off");
  RCP<EnhancedNumberValidator<int> > basicVali = 
  	rcp(new EnhancedNumberValidator<int>(1,10));
  boolValidatorDepList.set("do I have a validator?", 4, "does it have a validator?", basicVali);
  depSheet1->addDependency(RCP<BoolValidatorDependency>(
  	new BoolValidatorDependency(
	    My_deplist->getEntryRCP("Use Validator?"),
	    My_deplist->getEntryRCP("do I have a validator?"),
	    basicVali, 
	    RCP<ParameterEntryValidator>())));


  RCP<StringToIntegralParameterEntryValidator<int> >
    lowTempCheeseValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "PepperJack", "Swiss", "American" )
        ,"Cheese to Fondue"
        )
      );

  RCP<StringToIntegralParameterEntryValidator<int> >
    highTempCheeseValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "Munster", "Provalone", "Kurtis Super Awesome Cheese")
        ,"Cheese to Fondue"
        )
      );
  ParameterList&
    rangeValidatorDepList = My_deplist->sublist("Range Validator and NumberVisual Dependency List", false, "Range Validator and Number Visual Dependency testing list.");
  rangeValidatorDepList.set("Temperature",101.0, "The temperature of the fondue");
  rangeValidatorDepList.set("Cheese to Fondue", "Swiss", "The cheese we'll be using in our fondue pot.", lowTempCheeseValidator);
  RangeValidatorDependency<double>::RangeToValidatorMap tempranges;
  tempranges[std::pair<double,double>(100,200)] = lowTempCheeseValidator;
  tempranges[std::pair<double,double>(200,300)] = highTempCheeseValidator;
  RCP<RangeValidatorDependency<double> > cheeseTempDep = rcp(
  	new RangeValidatorDependency<double>(
	    My_deplist->getEntryRCP("Temperature"),
	    My_deplist->getEntryRCP("Cheese to Fondue"),
	    tempranges));
 
  depSheet1->addDependency(cheeseTempDep);
  
  RCP<SubtractionFunction<double> > fondueFunc = rcp(
    new SubtractionFunction<double>(100));
  RCP<NumberVisualDependency<double> > fondueDep = 
      RCP<NumberVisualDependency<double> >(new NumberVisualDependency<double>(
      My_deplist->getEntryRCP("Temperature"),
      My_deplist->getEntryRCP("Cheese to Fondue"),
      true,
      fondueFunc));
  depSheet1->addDependency(fondueDep);

  ParameterList&
    numberArrayLengthDepList = My_deplist->sublist("Number Array Length Dependency List", false, "Number Array Length ependecy testing list.");
  numberArrayLengthDepList.set("Array Length", 8, "array length setter");
  Array<double> variableLengthArray(10,23.0);
  RCP<EnhancedNumberValidator<double> > varLengthArrayVali = RCP<EnhancedNumberValidator<double> >(
  	new EnhancedNumberValidator<double>(10,50,4) );
  numberArrayLengthDepList.set("Variable Length Array", variableLengthArray, "variable length array",
  RCP<ArrayNumberValidator<double> >(new ArrayNumberValidator<double>(varLengthArrayVali)));

  RCP<NumberArrayLengthDependency<int, double> > arrayLengthDep(
    new NumberArrayLengthDependency<int, double>(My_deplist->getEntryRCP("Array Length"), 
    My_deplist->getEntryRCP("Variable Length Array")));
  depSheet1->addDependency(arrayLengthDep);



  ParameterList&
    numberValiAspDepList = My_deplist->sublist("Number Validator Aspect Dependency List", false, "Number Validator Aspect Dependency testing list.");
  RCP<EnhancedNumberValidator<int> > intVali2 = 
  	rcp(new EnhancedNumberValidator<int>(0,20));
  numberValiAspDepList.set("Int", 8, "Int tester", intVali2);
  numberValiAspDepList.set("Int2", 8, "int2 tester", intVali2);
  numberValiAspDepList.set("Int dependee", 1, "Int dependee");

  ParameterList&
    boolVisDepList = My_deplist->sublist("Bool Visual Dependency List", false, "Bool Visual Dependency testing list.");
  boolVisDepList.set("ShowPrecs", true, "Whether or not to should the Preciondtioner list");
  ParameterList&
    Prec_List0 = boolVisDepList.sublist("Preconditioner",false,"Sublist that defines the preconditioner.");
  Prec_List0.set("Type", "ILU", "The tpye of preconditioner to use");
  RCP<EnhancedNumberValidator<double> > droptolValidator = rcp(new EnhancedNumberValidator<double>(0,10,1e-3));
  Prec_List0.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.", droptolValidator);
  RCP<BoolVisualDependency> precDep1 = rcp(new BoolVisualDependency(
  My_deplist->getEntryRCP("ShowPrecs"),
  My_deplist->getEntryRCP("Preconditioner"),
  true));
  depSheet1->addDependency(precDep1);




ParameterList&
    stringVisDepList = My_deplist->sublist("String Visual Dependency List", false, "String Visual Dependency testing list.\nWorking June 29 2009");
  RCP<StringToIntegralParameterEntryValidator<int> >
    favCheeseValidator = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        tuple<std::string>( "Swiss", "American", "Cheder" )
        ,"Favorite Cheese"
        )
      );
   
   stringVisDepList.set(
   	"Favorite Cheese", "American", "Your favorite type of cheese", favCheeseValidator);
   RCP<EnhancedNumberValidator<int> > swissValidator = rcp(new EnhancedNumberValidator<int>(0,10));
   stringVisDepList.set("Swiss rating", 0, "How you rate swiss on a scale of 1 to 10", swissValidator);
   RCP<StringVisualDependency> swissDep1 = 
      RCP<StringVisualDependency>(new StringVisualDependency(
      My_deplist->getEntryRCP("Favorite Cheese"),
      My_deplist->getEntryRCP("Swiss rating"),
      "Swiss", 
      true));
   depSheet1->addDependency(swissDep1);





  RCP<SubtractionFunction<int> > intVisTester =
    rcp(new SubtractionFunction<int>(32));
  ParameterList&
    numberVisDepList = My_deplist->sublist(
      "Number Visual Dependency List", false, 
      "Number Visual Dependency testing list.");
  numberVisDepList.set("Ice", 50, "Ice stuff");
  numberVisDepList.set("Room Temp", 10, "Room temperature");
  RCP<NumberVisualDependency<int> > iceDep = 
      RCP<NumberVisualDependency<int> >(
        new NumberVisualDependency<int>(
      My_deplist->getEntryRCP("Room Temp"),
      My_deplist->getEntryRCP("Ice"),
      true,
      intVisTester));
  depSheet1->addDependency(iceDep);




  void (*myFunc)(RCP<const ParameterList>);
  myFunc = print;
  getInput(My_deplist, depSheet1, myFunc);

  std::cout << "Dep List: \n";
  writeParameterListToXmlOStream(*My_deplist, *out);

  std::cout << "Deps: \n";
  depSheet1->printDeps(std::cout);


}


}


int main(){
  Optika::RUN_OPTIKA_DEPENDENCY_AND_EXEC_TEST();
  return 0;
}


