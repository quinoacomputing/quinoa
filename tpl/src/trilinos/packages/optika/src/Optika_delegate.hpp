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
#ifndef OPTIKA_DELEGATE_HPP_
#define OPTIKA_DELEGATE_HPP_
#include <QItemDelegate>
#include <QModelIndex>
#include <QSize>
#include "Optika_ArrayWidget.hpp"

/*! \file Optika_delegate.hpp
    \brief The delegate used in the MVC
    framework for Optika.
*/

namespace Optika{

/**
 * \brief The delegate used for the Optika package. For non-documented functions please refer to the Qt API.
 */
class Delegate : public QItemDelegate{
	Q_OBJECT
public:
  /** \name Constructors */
  //@{

	/**
	 * \brief Constructs a Delegate.
	 * 
	 * @param parent The parent object of the Delegate.
	 */
	Delegate(QObject *parent = 0);

  //@}

  /** \name Overridden from QItemDelegate */
  //@{

  /** * \brief .  */
	QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
  /** * \brief .  */
	virtual void setEditorData(QWidget *editor, const QModelIndex &index) const;
  /** * \brief .  */
	void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;
  /** * \brief .  */
	virtual void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const;

  //@}

  /** \name Constant obtaining functions. */
  //@{

  /**
   * \brief Gets the value the delegate uses to represent 
   * "true" when constructing comboboxes for
   * boolean parameters.
   *
   * @return The value the delegate uses to represent 
   * "true" when constructing comboboxes for
   * boolean parameters.
   */
  static const QString& getBoolEditorTrue(){
    static const QString boolEditorTrue("true");
    return boolEditorTrue;
  }

  /**
   * \brief Gets the value the delegate uses to represent 
   * "false" when constructing comboboxes for
   * boolean parameters.
   *
   * @return The value the delegate uses to represent 
   * "false" when constructing comboboxes for
   * boolean parameters.
   */
  static const QString& getBoolEditorFalse(){
    static const QString boolEditorFalse("false");
    return boolEditorFalse;
  }

  //@}

private:
  /** \name Private Functions */
  //@{
  
  /**
   * \brief Creates an ArrayWidget to edit the parameter at the given QModelIndex.
   *
   * @param index Index of the parameter to edit.
   * @param type The template type of the array.
   * @param parent The parent widget.
   * @return An array widget with which to edit the parameter.
   */
	QWidget* getArrayEditor(const QModelIndex& index, QString type, QWidget *parent, bool isTwoD=false) const;

  /**
   * \brief Sets the data in an array widget to the current values found at index.
   *
   * @param editor ArrayWidget whose values are to be set.
   * @param type The template type of the array.
   * @param index The index of the parameter whose values are going to be
   * put in the ArrayWidget.
   */
  void setArrayWidgetData(QWidget* editor, QString type, const QModelIndex& index, bool isTwoD=false) const;

  /**
   * \brief Gets an array from an ArrayWidget and puts it into a QVariant.
   *
   * @param editor The ArrayWidget which has the array we want to put into
   * a QVariant.
   * @param type the template type of the array.
   * @return A Qvariant containing an array that is equal to the array in
   * the given ArrayWidget.
   */
  QVariant extractValueFromArray(QWidget* editor, QString type, bool isTwoD=false) const;

  //@}
};


}
#endif /* OPTIKA_DELEGATE_HPP_ */
