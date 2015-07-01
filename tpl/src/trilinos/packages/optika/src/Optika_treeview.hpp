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
#ifndef OPTIKA_TREEVIEW_HPP_
#define OPTIKA_TREEVIEW_HPP_
#include <QTreeView>
#include <QQueue>
#include "Optika_delegate.hpp"

/*! \file Optika_treeview.hpp
    \brief The view used in Optikas implementation
    of the MVC framework.
*/
namespace Optika{

class Delegate;
class TreeModel;

/**
 * \brief Class used to view TreeModels
 */
class TreeView : public QTreeView{
	Q_OBJECT
public:
  /** \name Public types */
  //@{

	/**
	 * \brief A pair representing an invalidIndex and why it's invalid
	 */
	typedef std::pair<QModelIndex, QString> invalidIndex;

  //@}

  /** \name Constructors */
  //@{

	/**
	 * \brief Constructs a TreeView.
	 * 
	 * @param treeModel The Tree Model being used with the TreeView.
	 * @param delegate The delegate to be used with the TreeView.
   * @param parent The parent widget.
	 */
	TreeView(TreeModel *treeModel, Delegate *delegate, QWidget* parent=0);

  //@}

public slots:

  /** \name Public Slots */
  //@{

	/**
	 * \brief Used to change the visiblity of a row from hidden to shown.
	 *
	 * @param row The row to be shwon.
	 * @param parent The parent of the item to be shown.
	 */
	void showRow(int row, const QModelIndex& parent);

	/**
	 * \brief Used to change the visiblity of a row from shown to hidden.
	 *
	 * @param row The row to be shwon.
	 * @param parent The parent of the item to be hidden.
	 */
	void hideRow(int row, const QModelIndex& parent);

	/**
	 * \brief Handles any badValue signals that might be emitted by the
	 * TreeModel.
	 *
	 * @param badValueIndex The index of the item with the bad value.
	 * @param A brief message explaining what happened to cause the
	 * treeitem to have an invalid value.
	 */
	void handleBadValue(QModelIndex badValueIndex, QString message);

	/**
	 * \brief Checks to see if there are any other invalid indicies.
	 * If there are, it dequeues the next invalidIndex from the 
	 * invalidIndicies queue and calls the handleBadValue function
	 * with it.
	 */
	void checkForOtherBadValues();

  //@}

private:
  /** \name Private Members */
  //@{
  
	/**
	 * \brief A Queue containing any invalid indicies that need to be
	 * delt with.
	 */
	QQueue<invalidIndex> invalidInicies; 

  //@}
};



}
#endif //OPTIKA_TREEVIEW_HPP_
