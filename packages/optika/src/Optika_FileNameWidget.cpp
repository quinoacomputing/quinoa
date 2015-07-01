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
#include "Optika_FileNameWidget.hpp"
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QFileDialog>


namespace Optika{


FileNameWidget::FileNameWidget(QString currentFileName, bool mustAlreadyExist, QWidget *parent)
	:QWidget(parent),
	currentFileName(currentFileName),
	mustAlreadyExist(mustAlreadyExist)
{
	QPushButton *changeButton = new QPushButton("Change Path",this);
	connect(changeButton, SIGNAL(clicked(bool)), this, SLOT(getNewFileName()));
	pathLabel = new QLabel(currentFileName,this);
	QVBoxLayout *layout = new QVBoxLayout(this);
	layout->addWidget(changeButton);
	layout->addWidget(pathLabel);
	setLayout(layout);
}

QString FileNameWidget::getCurrentFileName(){
	return currentFileName;
}

void FileNameWidget::setCurrentFileName(QString newName){
	currentFileName = newName;
	pathLabel->setText(newName);	
}

void FileNameWidget::getNewFileName(){
	QString defaultPath;
	if(currentFileName == ""){
		defaultPath = QDir::homePath();
	}
	else{
		defaultPath = currentFileName;
	}
	if(mustAlreadyExist){
    QString newFileName = QFileDialog::getOpenFileName(this, tr("File"), defaultPath);
    if(!newFileName.isNull()){
		  setCurrentFileName(newFileName);
    }
	}
	else{
		QString newFileName = QFileDialog::getSaveFileName(this, tr("File"), defaultPath);
    if(!newFileName.isNull()){
		  setCurrentFileName(newFileName);
    }
	}
}


}

