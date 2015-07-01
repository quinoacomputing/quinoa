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
#include "Optika_treeview.hpp"
#include <QMessageBox>
namespace Optika{


TreeView::TreeView(TreeModel *treeModel, Delegate *delegate, QWidget* parent):QTreeView(parent){
	setModel(treeModel);
	setItemDelegateForColumn(1, delegate);
	if(treeModel->hasDependencies()){
		connect(treeModel, SIGNAL(hideData(int, const QModelIndex&)), this, SLOT(hideRow(int, const QModelIndex&)));
		connect(treeModel, SIGNAL(showData(int, const QModelIndex&)), this, SLOT(showRow(int, const QModelIndex&)));
		connect(treeModel, SIGNAL(badValue(QModelIndex, QString)), this, SLOT(handleBadValue(QModelIndex, QString)));
		connect(delegate, SIGNAL(closeEditor(QWidget*, QAbstractItemDelegate::EndEditHint)), this, SLOT(checkForOtherBadValues()));
		treeModel->issueInitilizationSignals();
	}
	setAnimated(true);
	setAlternatingRowColors(true);
}

void TreeView::showRow(int row, const QModelIndex& parent){
	if(isRowHidden(row, parent)){
		setRowHidden(row, parent, false);
	}
}

void TreeView::hideRow(int row, const QModelIndex& parent){
	if(!isRowHidden(row, parent)){
		setRowHidden(row, parent, true);
	}
}

void TreeView::handleBadValue(QModelIndex badValueIndex, QString message){
	if(state() != EditingState && !isRowHidden(badValueIndex.row(), badValueIndex.parent())){
		QMessageBox::warning(this, "Bad parameter value", message);
		setCurrentIndex(badValueIndex);
		edit(badValueIndex);
	}
	else if(!isRowHidden(badValueIndex.row(), badValueIndex.parent())){
		invalidInicies.enqueue(invalidIndex(badValueIndex, message));
	}
}

void TreeView::checkForOtherBadValues(){
	if(invalidInicies.size() != 0){
		invalidIndex needsToBeEdited = invalidInicies.dequeue();
		handleBadValue(needsToBeEdited.first, needsToBeEdited.second);
	}
}

}

