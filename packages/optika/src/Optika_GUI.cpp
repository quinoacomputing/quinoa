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
#include <QApplication>
#include <QtGui>
#include <QString>
#include "Optika_GUI.hpp"
namespace Optika{

void getInput(RCP<ParameterList> validParameters, RCP<DependencySheet> dependencySheet, void (*customFunc)(RCP<const ParameterList>)){	
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow theWindow(validParameters, dependencySheet, customFunc);
		theWindow.show();
		a.exec();
	}
}

void getInput(
  const std::string& nameOfXmlFile,
  RCP<ParameterList>& userInput,
  void (*customFunc)(RCP<const ParameterList>))
{
  {
		using namespace Qt;
    RCP<DependencySheet> depSheet = rcp(new DependencySheet());
    userInput = getParametersFromXmlFile(nameOfXmlFile, depSheet);
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow theWindow(userInput, depSheet, customFunc);
		theWindow.show();
		a.exec();
	}
}

void getInputExtraOptions(
  RCP<ParameterList> validParameters,
  RCP<DependencySheet> dependencySheet,
  std::string styleSheetFilePath,
  std::string iconFilePath,
  void (*customFunc)(RCP<const ParameterList>)){	
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		std::string submitText = "Submit and save XML";
		std::string submitNoSaveText = "Submit";
		MetaWindow theWindow(validParameters, dependencySheet, customFunc, "", submitText, submitNoSaveText);
		if(iconFilePath != ""){
			QIcon windowIcon(QString::fromStdString(iconFilePath));
			a.setWindowIcon(windowIcon);
		}
		if(styleSheetFilePath != ""){
			QString str;
			QFile file(QString::fromStdString(styleSheetFilePath));
			if (file.open(QIODevice::ReadOnly | QIODevice::Text)){
				QTextStream in(&file);
				while (!in.atEnd()) {
					str += in.readLine();
				}
				a.setStyleSheet(str);
			}
		}

		theWindow.show();
		a.exec();
	}
}

void getInputExtraOptions(
  const std::string& nameOfXmlFile,
  RCP<ParameterList>& userInput,
  std::string styleSheetFilePath,
  std::string iconFilePath,
  void (*customFunc)(RCP<const ParameterList>))
{
  {
		using namespace Qt;
		RCP<DependencySheet> depSheet = rcp(new DependencySheet());
		userInput = getParametersFromXmlFile(nameOfXmlFile, depSheet);
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		std::string submitText = "Submit and save XML";
		std::string submitNoSaveText = "Submit";
		MetaWindow theWindow(userInput, depSheet, customFunc, "", submitText, submitNoSaveText);
		if(iconFilePath != ""){
			QIcon windowIcon(QString::fromStdString(iconFilePath));
			a.setWindowIcon(windowIcon);
		}
		if(styleSheetFilePath != ""){
			QString str;
			QFile file(QString::fromStdString(styleSheetFilePath));
			if (file.open(QIODevice::ReadOnly | QIODevice::Text)){
				QTextStream in(&file);
				while (!in.atEnd()) {
					str += in.readLine();
				}
				a.setStyleSheet(str);
			}
		}

		theWindow.show();
		a.exec();
	}
}


OptikaGUI::OptikaGUI(
  RCP<ParameterList> validParameters,
  RCP<DependencySheet> dependencySheet,
  void (*customFunc)(RCP<const ParameterList>)):
	validParameters(validParameters),
	dependencySheet(dependencySheet),
  customFunc(customFunc)
  {}

OptikaGUI::OptikaGUI(
  const std::string& xmlFileName,
  void (*customFunc)(RCP<const ParameterList>)):
  customFunc(customFunc)
{
  dependencySheet = rcp(new DependencySheet);
  validParameters = getParametersFromXmlFile(xmlFileName, dependencySheet);
}

void OptikaGUI::exec(){
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow *theWindow;
  	theWindow = new MetaWindow(validParameters, dependencySheet, customFunc);
		if(title != ""){
			theWindow->setWindowTitle(QString::fromStdString(title));
		}
		if(iconFilePath != ""){
			QIcon windowIcon(QString::fromStdString(iconFilePath));
			QApplication::setWindowIcon(windowIcon);
		}
		if(styleSheetFilePath != ""){
			QString str;
			QFile file(QString::fromStdString(styleSheetFilePath));
			if (file.open(QIODevice::ReadOnly | QIODevice::Text)){
				QTextStream in(&file);
				while (!in.atEnd()) {
					str += in.readLine();
				}
				a.setStyleSheet(str);
			}
		}
		if(aboutInfo != ""){
			theWindow->setAboutInfo(QString::fromStdString(aboutInfo));
		}
    theWindow->setActionButtonText(QString::fromStdString(actionButtonText));
		theWindow->show();
		theWindow->activateWindow();
		a.exec();
    delete theWindow;
	}
	
}

void OptikaGUI::setWindowTitle(const std::string& title){
	this->title = title;
}

void OptikaGUI::setWindowIcon(const std::string& filePath){
	this->iconFilePath = filePath;
}

void OptikaGUI::setStyleSheet(const std::string& filePath){
	this->styleSheetFilePath = filePath;
} 

void OptikaGUI::setCustomFunction(void (*customFunc)(RCP<const ParameterList>)){
	this->customFunc = customFunc;
}

void OptikaGUI::setAboutInfo(const std::string& aboutInfo){
	this->aboutInfo = aboutInfo;
}

void OptikaGUI::setActionButtonText(const std::string& text){
  actionButtonText = text;
}

std::string OptikaGUI::getWindowTitle(){
	return title;
}

std::string OptikaGUI::getWindowIcon(){
	return iconFilePath;
}

std::string OptikaGUI::getStyleSheet(){
	return styleSheetFilePath;
}

std::string  OptikaGUI::getAboutInfo(){
	return aboutInfo;
}

}

