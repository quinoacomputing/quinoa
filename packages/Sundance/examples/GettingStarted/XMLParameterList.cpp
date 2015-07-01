/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#include "Sundance.hpp"

int main(int argc, char** argv)
{
  try
  {
    /* Read the XML filename as a command-line option */
    string xmlFilename = "paramExample.xml";
    Sundance::setOption("xml-file", xmlFilename, "XML filename");
      
    /* Initialize */
    Sundance::init(&argc, &argv);

    /* Read a parameter list from the XML file */
    ParameterXMLFileReader reader(xmlFilename);
    ParameterList params = reader.getParameters();

    /* Get the parameters for the "Widget" sublist */
    const ParameterList& widget = params.sublist("Widget");
    Out::root() << "widget region label: " << widget.get<int>("Region") << endl;
    Out::root() << "widget material: " << widget.get<string>("Material") << endl;
    Out::root() << "widget density: " << widget.get<double>("Density") << endl;

    /* Get the parameters for the "Gizmo" sublist */
    const ParameterList& gizmo = params.sublist("Gizmo");
    Out::root() << "gizmo region label: " << gizmo.get<int>("Region") << endl;
    Out::root() << "gizmo material: " << gizmo.get<string>("Material") << endl;
    Out::root() << "gizmo density: " << gizmo.get<double>("Density") << endl;

    Sundance::passFailTest(true);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus();
}


