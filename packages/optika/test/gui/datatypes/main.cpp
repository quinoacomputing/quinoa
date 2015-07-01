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
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Optika_GUI.hpp"

namespace Optika{


void RUN_OPTIKA_DATA_TYPE_TESTS(){
  TEUCHOS_ADD_NUMBERTYPE_VALIDATOR_CONVERTERS(short)
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::FileNameValidator;
  using Teuchos::EnhancedNumberValidator;
  using Teuchos::StringValidator;
  using Teuchos::ArrayNumberValidator;
  using Teuchos::ArrayStringValidator;
  using Teuchos::ArrayFileNameValidator;
  using Teuchos::tuple;
  using Teuchos::rcp;

  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

  RCP<ParameterList> My_List2 = RCP<ParameterList>(new ParameterList);
  ParameterList&
    validatorList = My_List2->sublist("Validator List", false, "Validator testing\nWorking June 27th 2009");
  RCP<FileNameValidator> filnameVali = 
  	RCP<FileNameValidator>(new FileNameValidator);
  validatorList.set("filename", "", "filename tester", filnameVali);
  RCP<EnhancedNumberValidator<int> > intVali = 
  	rcp(new EnhancedNumberValidator<int>(0,10,2));
  validatorList.set("Int", 8, "Int tester", intVali);
  RCP<EnhancedNumberValidator<short> > shortVali = 
  	rcp(new EnhancedNumberValidator<short>(0,10,4));
  validatorList.set("Short", (short)4, "short tester", shortVali);
  RCP<EnhancedNumberValidator<float> > floatVali = 
  	rcp(new EnhancedNumberValidator<float>(0.0f,20.0f,1e-2f, 6));
  validatorList.set("Float", (float)4.5, "float tester", floatVali);
  RCP<EnhancedNumberValidator<double> > doubleVali = 
  	rcp(new EnhancedNumberValidator<double>(0,20,1e-2, 6));
  validatorList.set("Double", (double)4.5, "double tester", doubleVali);
  RCP<StringValidator> solverValidator2 = rcp(
      new StringValidator( tuple<std::string>( "GMRES", "CG", "TFQMR" )));
  validatorList.set(
    "Solver"
    ,"GMRES" // This will be validated by solverValidator right here!
    ,"The type of solver to use."
    ,solverValidator2
    );
  Array<std::string> validValues;
  validValues.append("value1");
  validValues.append("value2");
  validValues.append("value3");
  RCP<StringValidator> stringVali2 = RCP<StringValidator>(new StringValidator(validValues));
  validatorList.set("Easy String", "value1", "easy string validator tester", stringVali2);

  ParameterList&
    NoValiList = My_List2->sublist("No validator list",false,"sublist containing data types without validators on them for checking default behavior.");
  NoValiList.set("Int1", 8, "Int tester");
  NoValiList.set("Short1", (short)4, "short tester");
  NoValiList.set("Float1", (float)4.5, "float tester");
  NoValiList.set("Double1", (double)4.5, "double tester");
  NoValiList.set("Bool1", true);
  NoValiList.set("Bool", true);
  NoValiList.set("Free String", "fee");
  
  //Arrays
  RCP<StringValidator> easyStringValidator = RCP<StringValidator>(new StringValidator(tuple<std::string>("value1", "value2", "value3")));
  Array<int> intArray(10,0);
  Array<short> shortArray(10,3);
  Array<float> floatArray(10,4.4f);
  Array<double> doubleArray(10, 5.5);
  Array<std::string> stringArray(10,"CG");
  Array<std::string> easyStringArray(10, "value1");
  Array<std::string> freestringArray(10,"Blah");
  Array<std::string> filenameArray(3,"/net/home/f07/klnusbau/blah.txt");
  ParameterList&
  	ArrayList = My_List2->sublist("Arrays", false, "sublist containing arrays.");
  ArrayList.set("IntArray", intArray, "intarray tester", RCP<ArrayNumberValidator<int> >(new ArrayNumberValidator<int>(intVali)));
  ArrayList.set("ShortArray", shortArray, "shortarray tester", RCP<ArrayNumberValidator<short> >(new ArrayNumberValidator<short>(shortVali)));
  ArrayList.set("DoubleArray", doubleArray, "doublearray tester", RCP<ArrayNumberValidator<double> >(new ArrayNumberValidator<double>(doubleVali)));
  ArrayList.set("FloatArray", floatArray, "floatarray tester", RCP<ArrayNumberValidator<float> >(new ArrayNumberValidator<float>(floatVali)));
  ArrayList.set("StringArray", stringArray, "string tester", 
  RCP<ArrayStringValidator>(new ArrayStringValidator(solverValidator2))); 
  ArrayList.set("EasyStringArray", easyStringArray, "testing the easy validator", RCP<ArrayStringValidator>(new ArrayStringValidator(easyStringValidator)));
  ArrayList.set("FreeStringArray", freestringArray, "free string array tester");
  ArrayList.set("Filename Array", filenameArray, "filename array tester",
  	RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(filnameVali)));

  getInput(My_List2);
  writeParameterListToXmlOStream(*My_List2, *out);
  My_List2->print(
    std::cout,ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));



}

}

int main(){
  Optika::RUN_OPTIKA_DATA_TYPE_TESTS();
  return 0;
}
