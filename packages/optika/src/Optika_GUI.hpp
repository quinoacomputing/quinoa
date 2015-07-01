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
#ifndef OPTIKA_GUI_HPP_
#define OPTIKA_GUI_HPP_

#include "Optika_metawindow.hpp"

/*! \file Optika_GUI.hpp
    \brief A collection of functions and an
    Object that serve as the primary interface to the
    Optika package allowing application developers to
    obtain user input.*/

namespace Optika{

  //! @name Basic Input Getting Function
  //@{

  /**
   * \brief Retreives the input for a Teuchos Parameter List using a GUI. Note the 
   * Parameter List will be edited. All user input will be stored in it. 
   * Also runs the function specified whenever the user clicks the action
   * button and uses the specified dependency list.
   *
   * @param validParameters A list of parameters from which the users may 
   * specify values.
   * @param dependencySheet A sheet listing any dependencies between parameters
   * in the validParameters ParameterList.
   * @param customFunc Custom function to run whenever the user clicks the 
   * action button.
   */
  void getInput(
    RCP<ParameterList> validParameters,
    RCP<DependencySheet> dependencySheet=null,
    void (*customFunc)(RCP<const ParameterList>)=NULL);

  /**
   * \brief Reads in a set of parameters and dependencies from the specified xmlfile,
   * displays a GUI, and stores the users input in the sprecified
   * ParameterList. If a custom function is provided, it is run upon the user
   * clicking the action button.
   *
   * @param namOfXmlFile The name of the xml file from which parameters and
   * dependencies will be read in.
   * @param userInput A ParameterList into which all user input should be
   * stored.
   * @param customFunc A custom function for Optika to run upon the user 
   * clicking the action button.
   */
  void getInput(
    const std::string& nameOfXmlFile,
    RCP<ParameterList>& userInput,
    void (*customFunc)(RCP<const ParameterList>)=NULL);

  /**
   * \brief Retreives the input for a Teuchos Parameter List using a GUI. Note the 
   * Parameter List will be edited. All user input will be stored in it. 
   * Also runs the function specified whenever the user clicks the action
   * button and uses the specified dependency list.  One more button will be
   * available that will not ask user for saving the file.  A Qt stylesheet
   * file can be passed too.
   *
   * @param validParameters A list of parameters from which the users may 
   * specify values.
   * @param dependencySheet A sheet listing any dependencies between parameters
   * in the validParameters ParameterList.
   * @param styleSheetFilePath Path for Qt stylesheet
   * @param iconFilePath Path for application icon
   * @param customFunc Custom function to run whenever the user clicks the 
   * action button.
   */
  void getInputExtraOptions(
    RCP<ParameterList> validParameters,
    RCP<DependencySheet> dependencySheet=null,
    std::string styleSheetFilePath = "",
    std::string iconFilePath = "",
    void (*customFunc)(RCP<const ParameterList>)=NULL);

  /**
   * \brief Reads in a set of parameters and dependencies from the specified xmlfile,
   * displays a GUI, and stores the users input in the sprecified
   * ParameterList. If a custom function is provided, it is run upon the user
   * clicking the action button.  One more button will be available that will
   * not ask user for saving the file.  A Qt stylesheet file can be passed too.

   *
   * @param namOfXmlFile The name of the xml file from which parameters and
   * dependencies will be read in.
   * @param userInput A ParameterList into which all user input should be
   * stored.
   * @param styleSheetFilePath Path for Qt stylesheet
   * @param iconFilePath Path for application icon
   * @param customFunc A custom function for Optika to run upon the user 
   * clicking the action button.
   */
  void getInputExtraOptions(
    const std::string& nameOfXmlFile,
    RCP<ParameterList>& userInput,
    std::string styleSheetFilePath = "",
    std::string iconFilePath = "",
    void (*customFunc)(RCP<const ParameterList>)=NULL);


  //@}

/**
 * \brief A class that allows the user to create and customize their Optika GUI.
 */
class OptikaGUI{
public:
  /** \name Constructors */
  //@{

  /**
   * \brief Constructs an OptikaGUI object.
   *
   * @param validParameters A list of parameters from which the users may 
   * specify values.
   * @param dependencySheet A sheet listing any dependencies between parameters
   * in the validParameters ParameterList.
   * @param customFunc A custom function for Optika to run upon the user 
   * clicking the action button.
   */
  OptikaGUI(
    RCP<ParameterList> validParameters, 
    RCP<DependencySheet> dependencySheet=null,
    void (*customFunc)(RCP<const ParameterList>)=NULL);

  /**
   * \brief Constructs an OptikaGUI object.
   *
   * @param xmlFileName Name of an XML file describing the GUI.
   * @param customFunc A custom function for Optika to run upon the user 
   * clicking the action button.
   */
  OptikaGUI(const std::string& xmlFileName,
    void (*customFunc)(RCP<const ParameterList>)=NULL);

  //@}

  //! @name Execution Functions
  //@{
  
  /**
   * \brief Runs the GUI and gets the user input.
   */
  void exec();

  //@}

  //! @name Getters and Setters
  //@{
  
  /**
   * \brief Adds the information specified to the about dialog of the GUI.
   *
   * @param aboutInfo Information to be added to the about dialog of the GUI.
   */
  void setAboutInfo(const std::string& aboutInfo);

  /**
   * \brief Sets the text in the "action" button"
   *
   * @param text The text for the action button
   */
  void setActionButtonText(const std::string& text);

  /**
   * \brief Sets the title of the GUI window that is displayed to the user.
   *
   * @param title A string containing what the title of the GUI window 
   * should be.
   */
  void setWindowTitle(const std::string& title);

  /**
   * \brief Sets the window icon to the image specified in the filePath.
   *
   * @param filePath File path to the image that should be used as
   *  the window icon.
   */
  void setWindowIcon(const std::string& filePath);

  /**
   * \brief Sets the QT style sheet that should be used for the GUI.
   *
   * @param filePath File path to the QT style sheet to be used for
   * the GUI.
   */
  void setStyleSheet(const std::string& filePath);

  /**
   * \brief Sets the custom function to be used in the GUI. When ever the
   * user clicks the action button, this function will be run.
   *
   * @param The custom function to be run whenever the user clicks the action 
   * button.
   */
  void setCustomFunction(void (*customFunc)(RCP<const ParameterList>));

  /**
   * \brief Gets the window title.
   * 
   * @return A string containing the window title.
   */
  std::string getWindowTitle();

  /**
   * \brief Gets the file path describing the location of the file
   * being used for the window icon.
   *
   * @return The file path describing the location of the file
   * being used for the window icon.
   */
  std::string getWindowIcon();

  /**
   * \brief Gets the file path describing the location of the file
   * being used as the QT Style Sheet.
   *
   * @return The file path describing the location of the file
   * being used as the QT Style Sheet.
   */
  std::string getStyleSheet();

  /**
   * \brief Gets the information to be added to the about dialog of the GUI.
   *
   * @return the information to be added to the about dialog of the GUI.
   */
  std::string getAboutInfo();

  //@}

private:
  /** \name Private Members */
  //@{
  
  /**
   * \brief A list of parameters from which the users may specify values.
   */
  RCP<ParameterList> validParameters;

  /**
   * \brief A sheet listing any dependencies between parameters in the validParameters
   */
  RCP<DependencySheet> dependencySheet;

  /**
   * \brief A string containing the window title.
   */
  std::string title;

  /**
   * \brief File path to the image that should be used as the window icon.
   */
  std::string iconFilePath;

  /**
   * \brief File path to the QT style sheet to be used for the GUI.
   */
  std::string styleSheetFilePath;

  /**
   * \brief Information to be added to the about dialog of the GUI.
   */
  std::string aboutInfo;

  /**
   * \brief Text to display in the action button.
   */
  std::string actionButtonText;

  /**
   * \brief The custom function to be run whenever the user clicks the action button.
   */
  void (*customFunc)(RCP<const ParameterList>);

  //@}
};

}

#endif //OPTIKA_GUI_HPP_
