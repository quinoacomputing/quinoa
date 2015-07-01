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
#ifndef OPTIKA_ARRAYWIDGET_HPP_
#define OPTIKA_ARRAYWIDGET_HPP_

#include <QDialog>
#include <QModelIndex>
#include <QPushButton>
#include <QGridLayout>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QLineEdit>
#include <QGridLayout>
#include <QScrollArea>
#include <QLabel>
#include <vector>
#include "Optika_treemodel.hpp"
#include "Optika_FileNameWidget.hpp"
#include "Optika_ValidatorApplier.hpp"

/*! \file Optika_ArrayWidget.hpp
    \brief A collection of Widgets used to edit
    Arrays in  Optika.
*/

namespace Optika {

/**
 * \brief An Abstract base class for both 2D and 1D ArrayWidgets.
 *
 * Note the absence of the Q_OBJECT
 * macro. This is becuase classes using the Q_OBJECT macro can't be templated (bummer). The macro is therfore
 * present in the subclasses.
 */
template<class S>
class GenericArrayWidget : public QDialog{

public:

  /** \name Constructors */

  //@{

  /**
   * \brief Constructs a GenericArrayWidget.
   *
   * @param name The name of the parameter beting edited.
   * @param type The arrays template type.
   * @param validator The validator to be used on the Array.
   * @param parent The parent widget.
   */
	GenericArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0);
  //@}

	/**
	 * \brief Gets the type of array being edited.
	 *
	 * @return The type of array being edited.
	 */
	const QString getType() const {
		return type;
	}

  /**
   * \brief Returns the name of the parameter being edits.
   */
  const QString getName() const{
    return name;
  }

  /**
   * \brief Returns the validator being used on the array.
   */
  const RCP<const ParameterEntryValidator> getEntryValidator() const{
    return entryValidator;
  }

  /**
   * \brief called when the user is done entering data
   * into the widget. MUST BE IMPLEMENTED AS A SLOT IN
   * CONCRETE SUBCLASSES!
   */
  virtual void accept() =0;


protected: 

  /** @name Protected Functions */

  //@[

	/**
	 * \brief Sets up the layout for the arrayContainer, including adding what ever editing
	 * widget should be used for the particual type of array.
	 */
	virtual void setupArrayLayout(){
    if(arrayContainer->layout() == NULL){
      arrayContainer->setLayout(getArrayLayout());
    }
  }

  /**
   * \brief Get's the layout to be used for the array container in
   * the widget.
   *
   * @return The layout to be used for the array container in
   * the widget.
   */
  virtual QLayout* getArrayLayout() =0;  

	/**
   * \brief Gathers all the user inputed data and closes the dialog.
	 */
	virtual void doAcceptWork() =0;


  //@}

  /** @name Protected Members */
  //@{

	/**
	 * \brief The widget containing all of the editing widgets (e.g.
	 * QLineEdits, and QSpinBoxes) that comprise the array editor.
	 */
	QWidget *arrayContainer;

  //@}
	

private:

  /** @name Private Members */
  //@{

	/**
	 * \brief The type of array.
	 */
	QString type;

  /**
   * The name of the Parameter being edited.
   */
  QString name;

	/**
	 * \brief The validator being used on the array.
	 */
	RCP<const ParameterEntryValidator> entryValidator;	

  //@}
};

template<class S>
GenericArrayWidget<S>::GenericArrayWidget(
  QString name, 
  QString type, 
  const RCP<const ParameterEntryValidator> validator,
  QWidget *parent):
  QDialog(parent),
  type(type),
  name(name),
  entryValidator(validator)
{
	setModal(true);
	setSizeGripEnabled(true);
	arrayContainer = new QWidget(this);

	QScrollArea *scrollArea = new QScrollArea(this);
	scrollArea->setWidget(arrayContainer);
	scrollArea->setWidgetResizable(true);

	QPushButton *doneButton = new QPushButton(tr("Done"));
	QPushButton *cancelButton = new QPushButton(tr("Cancel"));
	connect(doneButton, SIGNAL(clicked(bool)), this, SLOT(accept()));
	connect(cancelButton, SIGNAL(clicked(bool)), this, SLOT(reject()));
	QGridLayout *layout = new QGridLayout(this);
	layout->addWidget(scrollArea,0,0,1,3);
	layout->addWidget(doneButton,1,2);
	layout->addWidget(cancelButton,1,1);

	this->setLayout(layout);

	setWindowTitle(name);
}

/**
 * \brief An abstract base class for all 2D Array Widets.
 *
 * \reference TwoDArray
 */
template <class S>
class Generic2DArrayWidget : public GenericArrayWidget<S>{

public:
  
  /** @name Constructors */
  //@{

  /**
   * \brief Constructs a Generic2DArrayWidget.
   *
   * @param name The name of the parameter beting edited.
   * @param type The arrays template type.
   * @param validator The validator to be used on the Array.
   * @param parent The parent widget.
   */
	Generic2DArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0);

  //@}

  /** @name Getters and Setters */
  //@{

  /**
   * \brief Initializes all of the data in the array widget
   * so when it pops up, the individual widgets are populate
   * with their current values in the array. Also, sets the
   * \c baseArray to \c array.
   *
   * @param array The array that should be used to populate
   * the individual widgets making up the ArrayWidget.
   */
  void initData(TwoDArray<S> array){
    baseArray = array;
    this->setupArrayLayout();
  }

  /**
   * \brief Retrieves all the data currently entered in each of the
   * individual widgets and compiles it into a TwoDArray.
   *
   * @return A TwoDArray representing the current values in the
   * idividual widgets.
   */
  TwoDArray<S> getArrayFromWidgets(){
    Teuchos::TwoDArray<QWidget*>::size_type numRows = 
      widgetArray.getNumRows()-1;
    Teuchos::TwoDArray<QWidget*>::size_type numCols = 
      widgetArray.getNumCols()-1;
    TwoDArray<S> toReturn(
      numRows, numCols);
    int numColsToIterate =0;
    for(int i=0; i<numRows; ++i){
      numColsToIterate = baseArray.isSymmetrical() ? 
        numCols-numRows+i : numCols;
      for(int j=0; j<numColsToIterate; ++j){
        toReturn(i,j) = getWidgetValue(i+1,j+1);
      }
    }
    toReturn.setSymmetrical(baseArray.isSymmetrical());
    return toReturn;
  }

  /**
   * \brief Gets the current value entred in the widget located at row,col.
   *
   * \return The value entered in the widget at row,col.
   */
  virtual S getWidgetValue(int row, int col) = 0;

  /**
   * \brief gets the data currently stored in baseArray.
   */
  inline TwoDArray<S> getData() const{
    return baseArray;
  }

  //@}

protected:
  
  /** @name Protected Functions */
  //@{

  /**
   * \brief Do all the things that need to be done when accept is called.
   *
   * Namely, set the base array to what is currently entered in the individual
   * widgets and call done.
   */
  void doAcceptWork(){
    baseArray.clear();
    baseArray = getArrayFromWidgets();
    this->done(QDialog::Accepted);
  }

  /**
   * Get a widget used for editing the value located at row,col
   * in the baseArray.
   *
   * @param row The row of the value for which an editor widget is desired.
   * @param col The column of the value for which and editor widget is desired.
   * @return An editor widget used to edit the value at row,co. in the
   * baseArray.
   */
  virtual QWidget* getEditorWidget(int row, int col) =0;

  //@}

  /** @name Protected members */
  //@{

  /** \brief An array containing the individual widgets which make up
   * the array widget editor. The widget at row,col is the widget which
   * will edit the value at row,col in the baseArray.
   */
  TwoDArray<QWidget*> widgetArray;

  /**
   * \brief The actual array data that the widget is editing. When the user
   * is finished editing, this TwoDArray is then populated with the values
   * they entered.
   */
  TwoDArray<S> baseArray;

  //@}

private:

  /** @name Private functions */
  //@{

  /**
   * \brief Retrieves the layout to be used in the array container.
   */
	QLayout* getArrayLayout();

};

template<class S>
Generic2DArrayWidget<S>::Generic2DArrayWidget(
  QString name, 
  QString type, 
  const RCP<const ParameterEntryValidator> validator,
  QWidget *parent):
  GenericArrayWidget<S>(name, type, validator, parent)
{}


template<class S>
QLayout* Generic2DArrayWidget<S>::getArrayLayout(){
 widgetArray = TwoDArray<QWidget*>(baseArray.getNumRows()+1, baseArray.getNumCols()+1);
 QGridLayout *widgetLayout = new QGridLayout;
  for(int i =0; i < baseArray.getNumCols(); ++i){
		widgetLayout->addWidget(new QLabel("Column: " +QString::number(i)),0,i+1,Qt::AlignLeft);
  }
  for(int i =0; i < baseArray.getNumRows(); ++i){
		widgetLayout->addWidget(new QLabel("Row: " +QString::number(i)),i+1,0,Qt::AlignLeft);
  }
  int numColsToIterate =0;
  for(int i =0; i < baseArray.getNumRows(); ++i){
    numColsToIterate = baseArray.isSymmetrical() ? 
      baseArray.getNumCols()-baseArray.getNumRows()+i : baseArray.getNumCols();
    for(int j =0; j < numColsToIterate; ++j){
		  QWidget* editorWidget = getEditorWidget(i,j);
		  widgetLayout->addWidget(editorWidget,i+1,j+1,Qt::AlignLeft);
		  widgetArray(i+1,j+1) = editorWidget;
    }
  }
  return widgetLayout;
}

/**
 * \brief A 2DArrayWidget used for editing arrays whose template type is int.
 */
class Int2DArrayWidget : public Generic2DArrayWidget<int>{
Q_OBJECT
public:

  /** @name Constructors */
  //@{

  /**
   * \brief Constructs an Int2DArrayWidget.
   *
   * @param name The name of the parameter beting edited.
   * @param type The arrays template type.
   * @param validator The validator to be used on the Array.
   * @param parent The parent widget.
   */
  Int2DArrayWidget(
    QString name,
    QString type,
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic2DArrayWidget<int>(name, type, validator, parent)
  {}


  /** @name Overridden from Generic2DArrayWidget */
  //@{

  /** \brief . */
  inline int getWidgetValue(int row, int col){
    return ((QSpinBox*)widgetArray(row,col))->value();
  }

protected:

  /** \brief . */
  QWidget* getEditorWidget(int row, int col){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<int> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const TwoDArrayValidator<EnhancedNumberValidator<int>, int> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<int>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray(row, col));
		return newSpin;
  }

  //@}

public slots:
  /** @name Overriden from GenericArrayWidget */
  //@{

  /** \brief . */
  void accept(){
    doAcceptWork();
  }

 //@} 

};

/**
 * \brief A 2DArrayWidget used for editing arrays whose template type is short.
 */
class Short2DArrayWidget : public Generic2DArrayWidget<short>{
Q_OBJECT
public:

  /** @name Constructors */
  //@{

  /**
   * \brief Constructs an Short2DArrayWidget.
   *
   * @param name The name of the parameter beting edited.
   * @param type The arrays template type.
   * @param validator The validator to be used on the Array.
   * @param parent The parent widget.
   */
  Short2DArrayWidget(
    QString name,
    QString type,
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic2DArrayWidget<short>(name, type, validator, parent)
  {}

  //@}

  /** @name Overridden from Generic2DArrayWidget */
  //@{

  /** \brief . */
  inline short getWidgetValue(int row, int col){
    return ((QSpinBox*)widgetArray(row,col))->value();
  }

protected:

  /** \brief . */
  QWidget* getEditorWidget(int row, int col){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<short> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const TwoDArrayValidator<EnhancedNumberValidator<short>, short> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<short>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray(row, col));
		return newSpin;
  }

  //@}

public slots:

  /** @name Overridden from GenericArrayWidget */
  //@{

  /** \brief . */
  void accept(){
    doAcceptWork();
  }

  //@}

};

/**
 * \brief A 2DArrayWidget used for editing arrays whose template type is double.
 */
class Double2DArrayWidget : public Generic2DArrayWidget<double>{
Q_OBJECT
public:

  /** @name Constructors */
  //@{

  /**
   * \brief Constructs an Double2DArrayWidget.
   *
   * @param name The name of the parameter beting edited.
   * @param type The arrays template type.
   * @param validator The validator to be used on the Array.
   * @param parent The parent widget.
   */
  Double2DArrayWidget(
    QString name,
    QString type,
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic2DArrayWidget<double>(name, type, validator, parent)
  {}

  //@}

  /** @name Overridden from Generic2DArrayWidget */
  //@{

  /** \brief . */
  inline double getWidgetValue(int row, int col){
    return ((QLineEdit*)widgetArray(row,col))->text().toDouble();
  }

protected:

  /** \brief . */
  QWidget* getEditorWidget(int row, int col){
		QLineEdit *newEdit = new QLineEdit(this);
		RCP<const EnhancedNumberValidator<double> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const TwoDArrayValidator<EnhancedNumberValidator<double>, double> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<double>::applyToLineEdit(validator, newEdit);
    newEdit->setText(QString::number(baseArray(row,col)));
		return newEdit;
  }

  //@}

public slots:

  /** @name Overridden from GenericArrayWidget */
  //@{

  void accept(){
    doAcceptWork();
  }

  //@}

};

/**
 * \brief A 2DArrayWidget used for editing arrays whose template type is float.
 */
class Float2DArrayWidget : public Generic2DArrayWidget<float>{
Q_OBJECT
public:

  /** @name Constructors */
  //@{

  /**
   * \brief Constructs an Float2DArrayWidget.
   *
   * @param name The name of the parameter beting edited.
   * @param type The arrays template type.
   * @param validator The validator to be used on the Array.
   * @param parent The parent widget.
   */
  Float2DArrayWidget(
    QString name,
    QString type,
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic2DArrayWidget<float>(name, type, validator, parent)
  {}

  //@}

  /** @name Overridden from Generic2DArrayWidget */
  //@{

  /** \brief . */
  inline float getWidgetValue(int row, int col){
    return ((QLineEdit*)widgetArray(row,col))->text().toDouble();
  }

protected:

  /** \brief . */
  QWidget* getEditorWidget(int row, int col){
		QLineEdit *newEdit = new QLineEdit(this);
		RCP<const EnhancedNumberValidator<float> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const TwoDArrayValidator<EnhancedNumberValidator<float>, float> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<float>::applyToLineEdit(validator, newEdit);
    newEdit->setText(QString::number(baseArray(row,col)));
		return newEdit;
  }

  //@}

public slots:

  /** @name Overridden from GenericArrayWidget */
  //@{

  void accept(){
    doAcceptWork();
  }

  //@}

};

/**
 * \brief A 2DArrayWidget used for editing arrays whose template type is std::string.
 */
class String2DArrayWidget : public Generic2DArrayWidget<std::string>{
Q_OBJECT
public:

  /** @name Constructors */
  //@{

  /**
   * \brief Constructs an String2DArrayWidget.
   *
   * @param name The name of the parameter beting edited.
   * @param type The arrays template type.
   * @param validator The validator to be used on the Array.
   * @param parent The parent widget.
   */
  String2DArrayWidget(
    QString name,
    QString type,
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic2DArrayWidget<std::string>(name, type, validator, parent)
  {}

  //@}

  /** @name Overridden from Generic2DArrayWidget */
  //@{

  /** \brief . */
  std::string getWidgetValue(int row, int col){
		if(is_null(getEntryValidator())){
       return ((QLineEdit*)widgetArray(row,col))->text().toStdString();
		}
		else if(!is_null(rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(getEntryValidator()))){
       return ((FileNameWidget*)widgetArray(row,col))->getCurrentFileName().toStdString();
		}
		else if(getEntryValidator()->validStringValues()->size() !=0){
       return  ((QComboBox*)widgetArray(row,col))->currentText().toStdString();
		}
		else{
       return  ((QLineEdit*)widgetArray(row,col))->text().toStdString();
		}
  }

protected:

  /** \brief . */
  QWidget* getEditorWidget(int row, int col){
    QString currentData = QString::fromStdString(baseArray(row,col));
		if(is_null(getEntryValidator())){
			return new QLineEdit(currentData,this);
		}
		else if(!is_null(rcp_dynamic_cast<const TwoDArrayValidator<FileNameValidator, std::string> >(getEntryValidator()))){
			return new FileNameWidget(
        currentData, 
        rcp_dynamic_cast<const TwoDArrayValidator<FileNameValidator, std::string> >(getEntryValidator())->getPrototype()->fileMustExist(), this);
		}
		else if(getEntryValidator()->validStringValues()->size() != 0){
			RCP<const Array<std::string> > options = getEntryValidator()->validStringValues();
			QComboBox *newCombo = new QComboBox(this);
			for(Array<std::string>::const_iterator itr = options->begin(); itr != options->end(); ++itr){
				newCombo->addItem(QString::fromStdString(*itr));
			}
      int selectedItem = newCombo->findText(currentData);
      newCombo->setCurrentIndex(selectedItem);
			return newCombo;
		}
		else{
			return new QLineEdit(currentData,this);
		}
  }

  //@}

public slots:
  /** @name Overridden from GenericArrayWidget */
  //@{

  /** \brief . */
  void accept(){
    doAcceptWork();
  }

  //@{

};

/**
 * \brief A templated abstract base class for all 1D array editing widgets. 
 *
 * \reference Array
 */ 
template <class S>
class Generic1DArrayWidget : public GenericArrayWidget<S>{
public:

  /** \name Constructor */
  //@{

  /**
   * \brief Constructs a Generic1DArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	Generic1DArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0);

  //@}
	
  //! @name Attribute/Query Methods 
  //@{


  /**
   * Return the array backing this widget.
   *
   * @param The array backing this widget.
   */
  const Array<S> getData() const{
    return baseArray;
  }

  //@}


  /** @name Miscellaneous */
  //@{

  /**
   * \brief Initializes all of the data in the array widget
   * so when it pops up, the individual widgets are populate
   * with their current values in the array. Also, sets the
   * \c baseArray to \c array.
   *
   * @param array The array that should be used to populate
   * the individual widgets making up the ArrayWidget.
   */
  void initData(Array<S> array){
    baseArray = array;
    this->setupArrayLayout();
  }

	/**
	 * \brief Gets the widget to be used as an editor for each entry in the array.
	 */
	virtual QWidget* getEditorWidget(int index) = 0;

  /**
   * \brief Get a new array reflecting the current values entered in the widgets.
   *
   * @return A new array reflecting the currecnt values entered in the widgets.
   */
  virtual Array<S> getArrayFromWidgets() = 0;


  //@}

protected:

  /** \name Protected types */
  //@{

	/**
	 * \brief Convienece typedef. Represents an array of QWidgets.
	 */
	typedef std::vector<QWidget*> WVector;

  //@}

  /** \name Protected members */
  //@{

	/**
	 * \brief Conatins the editing widgets (e.g. QLineEdits and QSpinBoxes) comprising the array editor.
	 */
	WVector widgetVector;

	/**
	 * \brief The array to be edited.
	 */
	Array<S> baseArray;
  
  //@}

  /** @name Overriden from GenericArrayWidget */
  //@{

  /** \brief . */
  void doAcceptWork();

  //@}

private:

  /** \name Private Functions */
  //@{

  QLayout* getArrayLayout();

  //@}
};

template<class S>
Generic1DArrayWidget<S>::Generic1DArrayWidget(
  QString name, 
  QString type, 
  const RCP<const ParameterEntryValidator> validator,
  QWidget *parent):
  GenericArrayWidget<S>(name, type, validator, parent)
{}


template<class S>
QLayout* Generic1DArrayWidget<S>::getArrayLayout(){
  QGridLayout *widgetLayout = new QGridLayout;
  for(int i=0; i<baseArray.size(); ++i){
	  widgetLayout->addWidget(new QLabel("Item: " +QString::number(i)),0,i,Qt::AlignLeft);
	  QWidget* editorWidget = getEditorWidget(i);
	  widgetLayout->addWidget(editorWidget,1,i,Qt::AlignLeft);
	  widgetVector.push_back(editorWidget);
  }
  return widgetLayout;
}

template<class S>
void Generic1DArrayWidget<S>::doAcceptWork(){
  baseArray.clear();
  baseArray = getArrayFromWidgets();
  this->done(QDialog::Accepted);
}

/**
 * \brief A 1D widget for editing Arrays of type int.
 */
class IntArrayWidget: public Generic1DArrayWidget<int>{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a IntArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	IntArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<int>(name, type, validator,parent){}

  //@}

  /** \name Overridden from Generic1DArrayWidget */
  //@{
  
  /**
   * \brief .
   */
	Array<int> getArrayFromWidgets(){
    Array<int> toReturn(widgetVector.size(), 0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QSpinBox*)widgetVector[i])->value();
    }
    return toReturn;
	}

private:

  /**
   * \brief .
   */
	QWidget* getEditorWidget(int index){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<int> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<int>, int> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<int>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

  //@}


public slots:

  /** \name Overridden from GenericArrayWidget */
  //@{

  /** * \brief .  */
  void accept(){
    doAcceptWork();
  }

  //@}

};

/**
 * \brief A widget for editing Arrays of type short.
 */
class ShortArrayWidget: public Generic1DArrayWidget<short>
{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a ShortArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	ShortArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<short>(name, type, validator,parent){}

  //@}

  /** @name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	Array<short> getArrayFromWidgets(){
    Array<short> toReturn(widgetVector.size(), 0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QSpinBox*)widgetVector[i])->value();
    }
    return toReturn;
	}


private:

  /** * \brief .  */
	QWidget* getEditorWidget(int index){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<short> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<short>, short> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<short>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

  //@}

public slots:
  /** \name Overridden from GenericArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}

  //@}
};

/**
 * \brief A widget for editing Arrays of type double.
 */
class DoubleArrayWidget: public Generic1DArrayWidget<double>
{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a DoubleArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	DoubleArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<double>(name, type, validator,parent){}

  //@}

  /** @name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	Array<double> getArrayFromWidgets(){
    Array<double> toReturn(widgetVector.size(), 0.0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QLineEdit*)widgetVector[i])->text().toDouble();
    }
    return toReturn;
  }

private:

  /** * \brief .  */
	QWidget* getEditorWidget(int index){
		QLineEdit *newEdit = new QLineEdit(this);
		RCP<const EnhancedNumberValidator<double> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<double>, double> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<double>::applyToLineEdit(validator, newEdit);
    newEdit->setText(QString::number(baseArray[index], 'g', ((QDoubleValidator*)newEdit->validator())->decimals()));
		return newEdit;
	}

  //@}
  
public slots:
  /** \name Overridden from GenericArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}

  //@}

};

/**
 * \brief A widget for editing Arrays of type short.
 */
class FloatArrayWidget: public Generic1DArrayWidget<float>
{
	Q_OBJECT
public:
  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a FloatArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	FloatArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<float>(name, type, validator,parent){}

  //@}
  
  /** @name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	Array<float> getArrayFromWidgets(){
    Array<float> toReturn(widgetVector.size(), 0.0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QLineEdit*)widgetVector[i])->text().toDouble();
    }
    return toReturn;
  }

private:

  /** * \brief .  */
	QWidget* getEditorWidget(int index){
		QLineEdit *newEdit = new QLineEdit(this);
		RCP<const EnhancedNumberValidator<float> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<float>, float> >(getEntryValidator(),true)->getPrototype();
		}
		ValidatorApplier<float>::applyToLineEdit(validator, newEdit);
    newEdit->setText(QString::number(baseArray[index], 'g', ((QDoubleValidator*)newEdit->validator())->decimals()));
		return newEdit;
	}


  //@}

public slots:

  /** \name Overridden from GenericArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}

  //@}

};

/**
 * \brief A widget for editing an array of strings
 */
class StringArrayWidget : public Generic1DArrayWidget<std::string>{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a StringArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
  StringArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<std::string>(name, type, validator,parent)
  {}

  //@}

  /** @name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	Array<std::string> getArrayFromWidgets(){
    Array<std::string> toReturn(widgetVector.size(), "");
    for(size_t i=0; i < widgetVector.size(); ++i){
			if(is_null(getEntryValidator())){
        toReturn[i] = ((QLineEdit*)widgetVector[i])->text().toStdString();
			}
			else if(!is_null(rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(getEntryValidator()))){
        toReturn[i] = ((FileNameWidget*)widgetVector[i])->getCurrentFileName().toStdString();
			}
			else if(getEntryValidator()->validStringValues()->size() !=0){
        toReturn[i] = ((QComboBox*)widgetVector[i])->currentText().toStdString();
			}
			else{
        toReturn[i] = ((QLineEdit*)widgetVector[i])->text().toStdString();
			}
		}
    return toReturn;
  }

private:

  /** * \brief .  */
	QWidget* getEditorWidget(int index){
    QString currentData = QString::fromStdString(baseArray[index]);
		if(is_null(getEntryValidator())){
			return new QLineEdit(currentData,this);
		}
		else if(!is_null(rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(getEntryValidator()))){
			return new FileNameWidget(
        currentData, 
        rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(getEntryValidator())->getPrototype()->fileMustExist(), this);
		}
		else if(getEntryValidator()->validStringValues()->size() != 0){
			RCP<const Array<std::string> > options = getEntryValidator()->validStringValues();
			QComboBox *newCombo = new QComboBox(this);
			for(Array<std::string>::const_iterator itr = options->begin(); itr != options->end(); ++itr){
				newCombo->addItem(QString::fromStdString(*itr));
			}
      int selectedItem = newCombo->findText(currentData);
      newCombo->setCurrentIndex(selectedItem);
			return newCombo;
		}
		else{
			return new QLineEdit(currentData,this);
		}
	}

  //@}

public slots:

  /** \name Overridden from GenericArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}

  //@}


};


} //end namespace
#endif //OPTIKA_ARRAYWIDGET_HPP_
