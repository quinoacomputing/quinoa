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
#include <Optika_GUI.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <QFileDialog>
#include <QApplication>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

  
static const std::string aboutInfo = 
"This is a generic ParameterList configurator for XML files "
"that conform to the standards of Teuchos' ParameterList serialization."
"\n\n"
"Developed by Kurtis Nusbaum and licensed under the LGPL\n\n";

const QString lastSaveDirTag("LastSaveDir");

QString getSettingsFileName(){
  QDir pldir = QDir::homePath() + "/.plconfigurator";
  if(!pldir.exists()){
    pldir.mkdir(pldir.absolutePath());
  }
  return QDir::homePath() + "/.plconfigurator/PLSettings.xml";
}


QString getLastSaveLocation(){
  QFile file(getSettingsFileName());
  QString toReturn = QDir::homePath();
  if(file.open(QIODevice::ReadOnly)){
    QXmlStreamReader xmlReader(&file);
		while(!xmlReader.isEndDocument()){
			if(xmlReader.isStartElement() && xmlReader.name() == lastSaveDirTag){
        toReturn = xmlReader.readElementText();
        break;
      }
			xmlReader.readNext();
    } 
    file.close();
  }
  return toReturn;
}

void saveLastSettings(QString lastSaveDir){
  QFile settings(getSettingsFileName());
  settings.open(QIODevice::WriteOnly);
  QXmlStreamWriter xmlWriter(&settings);
  xmlWriter.setAutoFormatting(true);
  xmlWriter.writeStartDocument();
  xmlWriter.writeStartElement("settings");
 	  xmlWriter.writeStartElement(lastSaveDirTag);
    xmlWriter.writeCharacters(lastSaveDir);
    xmlWriter.writeEndElement();
  xmlWriter.writeEndElement();
  xmlWriter.writeEndDocument();

}

void retreiveFileName(std::string& fileName){
		using namespace Qt;
    QString currentPath = QDir::currentPath();
		int argNum=1;
		char* args[1];
		std::string appName ="PLConfigurator";
		args[0] = &appName[0];
		QApplication a(argNum,args);
    QString caption("XML File to Configure");
    QString dirToStart = getLastSaveLocation();; 
    QString filter = "XML Files (*.xml)";
    QString retrieved = 
      QFileDialog::getOpenFileName(0, caption, dirToStart, filter);
    QFileInfo info(retrieved);
    saveLastSettings(info.dir().absolutePath());
    fileName = retrieved.toStdString();
}


int main(int argc, char* argv[]){
  using Teuchos::CommandLineProcessor;
  using Optika::OptikaGUI;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::DependencySheet;
  using Teuchos::getParametersFromXmlFile;

  bool verbose = true;
  bool success = true;
  try{
    std::string xmlFileName = "";
    CommandLineProcessor clp;
    clp.setOption(
      "xml-filename", 
      &xmlFileName,
      "The name of the xmlfile you wish to configure");
  
    if(xmlFileName == ""){
      retreiveFileName(xmlFileName);
      if(xmlFileName == ""){
        std::cout << "No file name provided. Exiting..." << std::endl;
        return 0;
      }
    }
  
    if(!QFile::exists(QString::fromStdString(xmlFileName))){
      std::cerr << "File does not exists!" << std::endl;
      std::cerr << "Bad file name: " << xmlFileName << "." << std::endl;
      return -1;
    }
  
  
    RCP<DependencySheet> deps = rcp(new DependencySheet);
  
    RCP<ParameterList> parameters = getParametersFromXmlFile(xmlFileName, deps);
  
    OptikaGUI og(parameters, deps);
    og.setAboutInfo(aboutInfo);
    og.setWindowTitle("PLEditor");
    og.setActionButtonText("Quit");
    og.setWindowIcon("plicon.svg");
    og.exec();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  return (success ? 0 : 1);

}
