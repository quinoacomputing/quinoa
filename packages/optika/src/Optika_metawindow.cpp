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
#include "Optika_metawindow.hpp"
#include "Optika_Version.hpp"
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox> 
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QtGui>
#include <QIcon>
#include <iostream>
#include <algorithm>
#include <QSettings>
namespace Optika{


const int numRecentDocuments = 7; 

SearchWidget::SearchWidget(TreeModel *treeModel, TreeView *treeView, QWidget *parent):
	QDialog(parent),
	treeModel(treeModel),
	treeView(treeView)
{
	matchesLabel = new QLabel(tr("Matches"));
	searchButton = new QPushButton(tr("Search"));
	connect(searchButton, SIGNAL(clicked(bool)), this, SLOT(search()));
	closeButton = new QPushButton(tr("Close"));
	connect(closeButton, SIGNAL(clicked(bool)), this, SLOT(hide()));
	nextButton = new QPushButton(tr("Next"));
	connect(nextButton, SIGNAL(clicked(bool)), this, SLOT(next()));
	previousButton = new QPushButton(tr("Previous"));
	connect(previousButton, SIGNAL(clicked(bool)), this, SLOT(previous()));
	searchTermsEdit = new QLineEdit(tr("Enter Search Terms Here"));
	QGridLayout *layout = new QGridLayout(this);
	layout->addWidget(searchTermsEdit,0,0);
	layout->addWidget(searchButton,0,1);
	layout->addWidget(nextButton,0,3);
	layout->addWidget(previousButton,0,2);
	layout->addWidget(closeButton,2,3);
	layout->addWidget(matchesLabel,3,0);
	setLayout(layout);
	nextButton->setDisabled(true);
	previousButton->setDisabled(true);
	setSizeGripEnabled(true);
	setWindowTitle(tr("Search..."));
}

void SearchWidget::search(){
	currentSearchResults = treeModel->match(treeModel->index(0,0,QModelIndex()), Qt::DisplayRole, searchTermsEdit->text(), 
	-1, Qt::MatchWrap | Qt::MatchContains | Qt::MatchRecursive);
	currentSearchResults = removeHiddenItems(currentSearchResults);
	currentSearchIterator = currentSearchResults.begin();	
	int searchSize = currentSearchResults.size();
	matchesLabel->setText("Matches ("+ QString::number(searchSize) + ")");
	if(searchSize <= 0){
		nextButton->setDisabled(true);
		previousButton->setDisabled(true);
	}
	else if(searchSize == 1){
		nextButton->setDisabled(true);
		previousButton->setDisabled(true);
		treeView->setCurrentIndex(*currentSearchIterator);
	}
	else{
		nextButton->setDisabled(false);
		previousButton->setDisabled(false);
		treeView->setCurrentIndex(*currentSearchIterator);
	}
}

void SearchWidget::next(){
	currentSearchIterator++;
	if(currentSearchIterator == currentSearchResults.end()){
		currentSearchIterator = currentSearchResults.begin();
	}
	treeView->setCurrentIndex(*currentSearchIterator);
}

void SearchWidget::previous(){
	currentSearchIterator--;
	if(currentSearchIterator == currentSearchResults.begin()-1){
		currentSearchIterator = currentSearchResults.end() -1;
	}
	treeView->setCurrentIndex(*currentSearchIterator);
}

QModelIndexList SearchWidget::removeHiddenItems(QModelIndexList& items){
	QModelIndexList toReturn;
	for(QModelIndexList::iterator it = items.begin(); it != items.end(); ++it){
		if(!treeView->isRowHidden(it->row(), it->parent())){
			toReturn.append(*it);
		}
	}
	return toReturn;
}

MetaWindow::MetaWindow(
  RCP<ParameterList> validParameters, 
  RCP<DependencySheet> dependencySheet, 
  void (*customFunc)(RCP<const ParameterList>), 
  QString fileName,
  const std::string actionButtonText,
  const std::string actionNoSaveButtonText)
{
	model = new TreeModel(validParameters, dependencySheet, fileName, this);
	initilization(customFunc, actionButtonText, actionNoSaveButtonText);
} 


MetaWindow::~MetaWindow(){
	saveSettings();
}

void MetaWindow::closeEvent(QCloseEvent *event){
	if(!model->isSaved()){
		if(saveCurrentUnsavedFile()){
			saveSettings();
			event->accept();
		}
		else{
			event->ignore();
		}
	}
	else{
		saveSettings();
		event->accept();
	}
}

void MetaWindow::initilization(
  void (*customFunc)(RCP<const ParameterList>),
  const std::string actionButtonText,
  const std::string actionNoSaveButtonText)
{
	this->customFunc = customFunc;
	delegate = new Delegate(this);
	view = new TreeView(model, delegate, this);
	view->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::SelectedClicked);
	searchWidget = new SearchWidget(model, view, this);
	searchWidget->hide();
  actionButton = NULL;
  actionNoSaveButton = NULL;
  if(actionButtonText != ""){
    actionButton = 
      new QPushButton(QString::fromStdString(actionButtonText), this); 
  }
  else{
	  actionButton = new QPushButton(tr("Submit"), this);
  }
  if(actionNoSaveButtonText != ""){
      actionNoSaveButton = new QPushButton(QString::fromStdString(actionNoSaveButtonText), this); 
  }
	QWidget *centerWidget = new QWidget(this);
	QGridLayout *centerWidgetLayout = new QGridLayout(centerWidget);
	centerWidgetLayout->addWidget(view,0,0);
	connect(actionButton, SIGNAL(clicked(bool)), this, SLOT(doAction()));
	if(actionNoSaveButton){
	    centerWidgetLayout->addWidget(actionButton,1,0,Qt::AlignLeft);
		connect(actionNoSaveButton, SIGNAL(clicked(bool)), this, SLOT(doActionNoSave()));
		centerWidgetLayout->addWidget(actionNoSaveButton,1,0,Qt::AlignRight);
	}
    else
	  centerWidgetLayout->addWidget(actionButton,1,0,Qt::AlignRight);
	centerWidget->setLayout(centerWidgetLayout);
	setCentralWidget(centerWidget);

	createActions();
	createMenus();
	resize(800,600);
	currentLoadDir = QDir::homePath();
	currentSaveDir = QDir::homePath();
	loadLastSettings();
	setWindowTitle(tr("Parameter Input"));
	view->show();
	view->header()->resizeSections(QHeaderView::ResizeToContents);
	view->header()->setMovable(false);
}

void MetaWindow::createMenus(){
//	recentMenu = new QMenu(tr("Recent Solvers"));
//	QAction *noRecentAct = new QAction(tr("No Recent Documents"),this);
//	noRecentAct->setEnabled(false);
//	recentMenu->addAction(noRecentAct);
	fileMenu = menuBar()->addMenu(tr("File"));
	fileMenu->addAction(resetAct);
	//fileMenu->addMenu(recentMenu);
	fileMenu->addSeparator();
	fileMenu->addAction(saveAct);
	fileMenu->addAction(saveAsAct);
	fileMenu->addAction(loadAct);
	fileMenu->addSeparator();
	fileMenu->addAction(quitAct);
	helpMenu = menuBar()->addMenu(tr("Help"));
	helpMenu->addAction(aboutAct);
	helpMenu->addAction(searchAct);
}

void MetaWindow::createActions(){
	resetAct = new QAction(tr("&Reset"),this);
	resetAct->setShortcut(tr("Ctrl+R"));
	resetAct->setStatusTip(tr("Reset the list to its original state."));
	connect(resetAct, SIGNAL(triggered()), this, SLOT(resetModel()));

	saveAct = new QAction(tr("&Save"),this);
	saveAct->setShortcut(tr("Ctrl+S"));
	saveAct->setStatusTip(tr("Save the current file."));
	connect(saveAct, SIGNAL(triggered()), this, SLOT(saveFile()));

	saveAsAct = new QAction(tr("Save As..."),this);
	saveAsAct->setStatusTip(tr("Save the current file to a specified file name."));
	connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveFileAs()));

	loadAct = new QAction(tr("&Load"),this);
	loadAct->setShortcut(tr("Ctrl+L"));
	loadAct->setStatusTip(tr("Load input file"));
	connect(loadAct, SIGNAL(triggered()), this, SLOT(loadFile()));

	quitAct = new QAction(tr("&Quit"),this);
	quitAct->setShortcut(tr("Ctrl+Q"));
	quitAct->setStatusTip(tr("Quit"));
	connect(quitAct, SIGNAL(triggered()), this, SLOT(close()));

	aboutAct = new QAction(tr("About"),this);
	searchAct = new QAction(tr("Search For Parameter/Parameter List"), this);
	searchAct->setToolTip("Search for a particular Parameter or ParameterList");
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(showAbout()));
	connect(searchAct, SIGNAL(triggered()), this, SLOT(initiateSearch()));
}

void MetaWindow::load(){
	QString fileName = QFileDialog::getOpenFileName(this, tr("Load..."), currentLoadDir, tr("Xml (*.xml)"));
	if(fileName != ""){
		model->readInput(fileName);
		currentLoadDir = fileName.section("/",0,-2);
		addRecentDocument(fileName);
	}
}

void MetaWindow::loadLastSettings(){
  QSettings settings(QSettings::UserScope, "Sandia", "Optika");
  currentSaveDir = 
    settings.value(lastSaveDirSetting(), currentSaveDir).toString();
  currentLoadDir = 
    settings.value(lastLoadDirSetting(), currentLoadDir).toString();
  resize(settings.value(xresSetting(), width()).toInt(),
    settings.value(yresSetting(), height()).toInt());
  move(settings.value(xposSetting(), x()).toInt(),
    settings.value(yposSetting(), y()).toInt());
}


void MetaWindow::saveSettings(){
  QSettings settings(QSettings::UserScope, "Sandia", "Optika");

  settings.setValue(lastSaveDirSetting(), currentSaveDir);
  settings.setValue(lastLoadDirSetting(), currentLoadDir);
  settings.setValue(xresSetting(), width());
  settings.setValue(yresSetting(), height());
  settings.setValue(xposSetting(), x());
  settings.setValue(yposSetting(), y());
}
	

void MetaWindow::addRecentDocument(QString recentDocument){
	recentDocsList.prepend(recentDocument);
	if(recentDocsList.size() > numRecentDocuments){
		recentDocsList.removeLast();
	}
//	updateRecentDocsMenu();
}

void MetaWindow::updateRecentDocsMenu(){
	recentMenu->clear();
	for(int i=0; i<recentDocsList.size(); ++i){
		QAction *recentDocAct = new QAction(recentDocsList.at(i).section("/",-1,-1),this);
		connect(recentDocAct, SIGNAL(triggered()), this, SLOT(loadRecentDoc()));
		recentMenu->addAction(recentDocAct);
	}
}

void MetaWindow::resetModel(){
	if(!model->isSaved()){
		saveCurrentUnsavedFile();
	}
	model->reset();
}

bool MetaWindow::saveFileAs(){
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save To..."), currentSaveDir, tr("XML (*.xml)"));
	if(fileName.toStdString() != ""){
		if(!fileName.endsWith(".xml")){
			fileName = fileName.append(".xml");
		}
		if(model->writeOutput(fileName)){
			currentSaveDir = fileName.section("/",0,-2);
			addRecentDocument(fileName);		
			return true;
		}
	}
	return false;
}

void MetaWindow::saveFile(){
	QString currentFileName = model->getSaveFileName();
	if(currentFileName != ""){
		model->writeOutput(currentFileName);
	}
	else{
		saveFileAs();
	}
}


void MetaWindow::loadFile(){
	if(!model->isSaved()){
		saveCurrentUnsavedFile();
	}
	load();
}

bool MetaWindow::saveCurrentUnsavedFile(){
		QMessageBox saveQuestion(QMessageBox::Question, 
		tr("Save?"),
		tr("These choices have not been saved since you last made changes. Would you like to save them now?"), 
		QMessageBox::Yes | QMessageBox::No,
		this);
		saveQuestion.setDefaultButton(QMessageBox::Yes);
		int shouldSave = saveQuestion.exec(); 
		if(shouldSave == QMessageBox::Yes){
			return saveFileAs();
		}
		return true;
}

void MetaWindow::loadRecentDoc(){
	QString docName = dynamic_cast<QAction*>(sender())->text();
	int i =0;
	for(; i<recentDocsList.size();++i){
		if(recentDocsList.at(i).contains(docName)){
			break;
		}
	}
	if(!model->isSaved()){
		if(saveCurrentUnsavedFile()){
			model->readInput(recentDocsList.at(i));
		}
	}
}

void MetaWindow::showAbout(){
	QString aboutString = aboutInfo + "\nThis input obtainer was generated by Kurtis Nusbaum's Optika package, part of the Trilinos Project.\n\nVersion: " + QString::fromStdString(Optika_Version()) + "\nWebsite: trilinos.sandia.gov/packages/optika\nLicense: LGPL\nContact: klnusbaum@gmail.com";
	QMessageBox::about(this,
	"Optika Input Obtainer\n",
	aboutString);
}

void MetaWindow::initiateSearch(){
	searchWidget->show();
}

void MetaWindow::doAction(){
	if(customFunc == 0){
		close();
	}
	else{
		(*customFunc)(model->getCurrentParameters());	
	}
}

void MetaWindow::doActionNoSave(){
	model->setIsSaved();
    doAction();
}

void MetaWindow::setAboutInfo(QString aboutInfo){
	this->aboutInfo = aboutInfo;
}

QString MetaWindow::getAboutInfo(){
	return aboutInfo;
}

void MetaWindow::setActionButtonText(QString newText){
  actionButton->setText(newText);
}

QString MetaWindow::getActionButtonText(){
  return actionButton->text();
}


} //namespace Optika

