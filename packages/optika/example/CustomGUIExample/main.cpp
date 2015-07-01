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

int main(int argc, char* argv[])
{
 /*
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  * !!!!!!!!!!!!!!!!              ATTENTION              !!!!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!   PLEASE VIEW THE BASIC EXAMPLE FIRST BEFORE READING THIS       !!!!!!
  * !!!!   EXAMPLE. IT PROVIDES FUNDAMENTAL KNOWLEDGE THAT WILL BE VERY  !!!!!! 
  * !!!!   HELPFUL IN UNDERSTANDING THIS EXAMPLE.                        !!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  */ 

  /*
   * Defautls not good enought for you, eh? Well I guess That's understandable 
   * and why we've given you a few ways to customize your GUI as you see fit. 
   * The purpose of this example is to demonstrate those capabilities.
   */

  /* 
   * So let's setup a list of parameters to obtain. They could be anything 
   * really. We'll just reuse the list we created in the Basic Example 
   * (which you've read already because you're smart and headed the gigantic
   * warning I put at the top of this example.
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
    "The tolerance below which entries from the"
    "factorization are left out of the factors.");

  /**
   * Now here's were things get a little differnent. Instead of just calling 
   * Optika::getInput(), we're actually going to create and OptikaGUI object.
   * It will be the vehical through which we customize the GUI. We'll pass it
   * the ParameterList in which we want it to story user input.
   */
   Optika::OptikaGUI myGUI(My_List);

  /**
   * Now we can start configuring what our GUI will look like. Let's start 
   * with the window title.
   */
   myGUI.setWindowTitle("My Custom GUI");

  /**
   * We can set the information that will be displayed in the about dialog for 
   * the gui.
   */
   myGUI.setAboutInfo("This is a little GUI I made.");

  /**
   * If you have an icon you'd like to use as the WindowIcon you can do that 
   * too. Just specify the path to the file containig the icon. Supported 
   * formats will vary from system to system and your QT installation, but the 
   * following should always work:
   *	-BMP   -GIF  -JPG  -JPEG
	 *  -MNG   -PNG  -PBM  -PGM
	 *  -PPM   -TIFF -XBM  -XPM
	 *  -SVG
   */
   myGUI.setWindowIcon("myIcon.png");

  /**
   * Now if you really wanna dive into your GUI design you can use
   * QT style sheets. Since optika is build on top of QT, you can
   * use QT style sheets to style the various widgets used by
   * Optika. The main widges Optika uses thay you'll probably 
   * wanna style are:
      -QTreeView -QDialog -QMainWindow
      -QPushButton -QSpinBox -QMenu
	    -QMenuBar
   * You might need to look at some of the Optika source code to really
   * get fine-grained styling control. Once your stylesheet is made,
   * all you need to do is specify the filepath to the Optika_GUI
   * object. Also note the style sheet provided in this example is 
   * exceptionally ugly.
   */
   myGUI.setStyleSheet("myStyleSheet.qss");

  /**
   * Now that we're all ready to go, we just call the exec funtion on
   * our Optika_GUI object. This will get the user input and put it in our
   * pramater list.
   */
   myGUI.exec();

   /**
    * That's it! You can make even more awesome GUIs now!.
	  */
  return 0;
}

