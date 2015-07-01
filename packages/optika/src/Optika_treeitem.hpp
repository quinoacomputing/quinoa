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
#ifndef OPTIKA_TREEITEM_HPP_
#define OPTIKA_TREEITEM_HPP_

#include <QList>
#include <QVariant>
#include "Optika_ArrayHelperFunctions.hpp"

/*! \file Optika_treeitem.hpp
    \brief A item in the treemodel.
*/

namespace Optika{
/**
 * \brief The TreeItem class is the item class used by the TreeModel class.
 */
class TreeItem{
public:
  /** \name Constructors/Destructor */
  //@{

	/**
	 * \brief Constructs a TreeItem object.
	 *
   * @param name The name of the parameter entry.
	 * @param parameterEntry The ParameterEntry this TreeItem is ment to represent.
	 * @param parent The parent TreeItem.
   * @param isHeader Whether or not his treeitem represents a "header" tree item.
	 */
	TreeItem(const QString& name, RCP<ParameterEntry> parameterEntry, TreeItem *parent = 0, bool isHeader=false);

	/**
	 * \brief Deconstrcutor for the TreeItem.
	 */
	~TreeItem();

  //@}

  //! @name Debugging fucntions
  //@{

	/**
	 * \brief Prints out the values in the TreeItem.
	 */
	void printOut() const;

  //@}

  //! @name Getters and Setters
  //@{
  
	/**
	 * \brief Returns the child treeitem in the row specified by the row argument.
	 *
	 * @param row The row in which the child is in.
	 * @return The child TreeItem located in the row.
	 */
	TreeItem *child(int row);

  /**
   * \brief Gets the ParameterEntry associated with this TreeItem.
   *
   * @return The ParameterEntry associated with this TreeItem. If this
   * tree item does not have a parameterEntry, null is returned.
   */
  inline RCP<const ParameterEntry> getEntry() const{
    return parameterEntry.getConst();
  }

  /**
   * \brief Returns whether or not this TreeItem has a ParameterEntry associated
   * with it.
   */
  inline bool hasEntry() const{
    return parameterEntry != null;
  }

	/**
	 * \brief Gets the number of child nodes this item has.
	 *
	 * @return The number of child nodes this item has.
	 */
	int childCount() const;

	/**
	 * \brief Gets a list of all the child items.
	 *
	 * @return A list of all child items.
	 */
	const QList<TreeItem*> getChildItems();

	/**
	 * \brief How man columns the TreeItem has. Should always be 3.
	 *
	 * @return The number of columns the TreeItem has.
	 */
	int columnCount() const;

	/**
	 * \brief Returns the data located in a particular column.
	 *
	 * @param column The column of the desired data.
	 * @param role The role of the data.
	 * @return The data located in a particular column.
	 */
	QVariant data(int column, int role = Qt::DisplayRole) const;

	/**
	 * \brief Gets the parent TreeItem
	 *
	 * @return The parent TreeItem.
	 */
	TreeItem *parent();

	/**
	 * \brief Returns the row in which this TreeItem is located.
	 * 
	 * @return The row in which this TreeItem is located.
	 */
	int row() const;

	/**
	 * \brief Determines whether or not the current value associated with the
	 * TreeItem is valid.
	 *
	 * @return True if the value is valid, false otherwise.
	 */
	bool hasValidValue() const;

  /**
   * \brief Gets a message desribing the error with the current value.
   *
   * This funciton returns the error message generated by the validator
   * should this treeitem's current value be invalid. If the current value
   * is actually valid then this function simply returns an empty string.
   * If there is no validator on this treeitem then an empty string is
   * returned.
   *
   * @return The error message generated by the validator. If the current
   * value is valid or there is no validator then the string is empty.
   */
  QString getCurrentInvalidValueMessage() const;

	/**
	 * \brief Appends a child TreeItem to the TreeItem
	 * 
	 * @param child The child item to be appended.
	 */
	void appendChild(TreeItem *child);

	/**
	 * \brief Changes the value of the TreeItem. Should only be used with TreeItems that represent Parameters.
	 *
	 * @param value The new value to be assigned to the TreeItem.
	 */
	bool changeValue(QVariant value);

	/**
	 * \brief Sets the validator for the parameter the TreeItem represents.
	 *
	 * @param validator The validator which the parameter should be given.
	 */
	void setValidator(RCP<const ParameterEntryValidator> validator);

  /**
   * \brief Gets the type id to be used for the TreeItem.
   *
   * @param param The parameter who's type is in question.
   * @return The type id for the parameter.
   */
  static QString getTypeId(const RCP<const ParameterEntry> parameter);


  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief the name of the tree item.
   */
  const QString name;

  /**
   * \brief The type id associated with this TreeItem.
   */
  QString myTypeId;

	/**
	 * \brief The childitems of the TreeItem.
	 */
	QList<TreeItem*> childItems;

	/**
	 * \brief The parent TreeItem.
	 */
	TreeItem *parentItem;

	/**
	 * \brief The ParameterEntry being represented by the TreeItem.
	 */
	RCP<ParameterEntry> parameterEntry;

	/**
	 * \brief The docString for the TreeItem.
	 */
	QString docString;

  /**
   * \brief Whether or not this is a header treeitem.
   */
  bool isHeader;

  //@}

  /** \name Private Functions */
  //@{
  
	/**
	 * \brief Changes the value of an array.
	 *
	 * @param value A string representing the value of the array.
	 * @param type The type of array.
	 */
	void changeValueForArray(QVariant value, QString type, bool twoD=false);

  //@}

};



}

#endif /* OPTIKA_TREEITEM_HPP_ */
