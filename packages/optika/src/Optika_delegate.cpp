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
#include <QtGui>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QComboBox>
#include <QFileDialog>
#include "Optika_delegate.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace Optika{

Delegate::Delegate(QObject *parent):QItemDelegate(parent){}

QWidget* Delegate::createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/ , const QModelIndex &index ) const{
	QWidget *editor = 0;
	if(index.column() != 1)
		return editor;

	RCP<const ParameterEntryValidator> paramValidator = ((TreeModel*)(index.model()))->getValidator(index);
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);

	if(itemType == intId){
		editor = new QSpinBox(parent);
		RCP<const EnhancedNumberValidator <int> > intValidator;
		if(!is_null(paramValidator)){
			intValidator = rcp_dynamic_cast<const EnhancedNumberValidator<int> >(paramValidator);
		}
		ValidatorApplier<int>::applyToSpinBox(intValidator, (QSpinBox*)editor);
	}
	else if(itemType == shortId){
		editor = new QSpinBox(parent);
		RCP<const EnhancedNumberValidator<short> > shortValidator;
		if(!is_null(paramValidator)){
			shortValidator = rcp_dynamic_cast<const EnhancedNumberValidator<short> >(paramValidator);
		}
		ValidatorApplier<short>::applyToSpinBox(shortValidator, (QSpinBox*)editor);
	}
/*	else if(itemType == longlongId){
		editor = new QwwLongSpinBox(parent);
		RCP<const EnhancedNumberValidator<long long> > longlongValidator;
		if(!is_null(paramValidator)){
			longlongValidator = rcp_dynamic_cast<const EnhancedNumberValidator<long long> >(paramValidator);
		}
		EnhancedNumberValidator<long long>::applyToSpinBox(longlongValidator, (QDoubleSpinBox*)editor);
	}*/
	else if(itemType == doubleId){
		//editor = new QLineEdit(parent);
		editor = new QLineEdit(parent);
		RCP<const EnhancedNumberValidator<double> > doubleValidator;
		if(!is_null(paramValidator)){
			doubleValidator = rcp_dynamic_cast<const EnhancedNumberValidator<double> >(paramValidator);
		}
		ValidatorApplier<double>::applyToLineEdit(doubleValidator, (QLineEdit*)editor);
	}
	else if(itemType == floatId){
		editor = new QLineEdit(parent);
		RCP<const EnhancedNumberValidator<float> > floatValidator; 
		if(!is_null(paramValidator)){
			floatValidator = rcp_dynamic_cast<const EnhancedNumberValidator<float> >(paramValidator);
		}
		ValidatorApplier<float>::applyToLineEdit(floatValidator, (QLineEdit*)editor);
	}
	else if(itemType == boolId){
		editor = new QComboBox(parent);
		static_cast<QComboBox*>(editor)->addItem(getBoolEditorTrue());
		static_cast<QComboBox*>(editor)->addItem(getBoolEditorFalse());
	}
	else if(itemType == stringId){
		if(is_null(paramValidator)){
			editor = new QLineEdit(parent);
		}
		else if(!is_null(rcp_dynamic_cast<const FileNameValidator>(paramValidator))){
			QString paramName = 
				((TreeModel*)(index.model()))->data(index.sibling(index.row(), 0),Qt::DisplayRole).toString();
			QString currentPath = ((TreeModel*)(index.model()))->data(index,Qt::DisplayRole).toString();
			if(currentPath.size() == 0){
				currentPath = QDir::homePath();
			}
			// Hack.
			if(paramName.indexOf("Input Directory", 0, Qt::CaseInsensitive) > -1 ||
			   paramName.indexOf("Output Directory", 0, Qt::CaseInsensitive) > -1){
				QString dirname = QFileDialog::getExistingDirectory(parent, paramName,
					currentPath, QFileDialog::ShowDirsOnly);
				if(dirname != ""){
					((TreeModel*)(index.model()))->setData(index, dirname);
				}
			}
			else{
				QString filename;
				if(rcp_dynamic_cast<const FileNameValidator>(paramValidator)->fileMustExist()){
					filename = QFileDialog::getOpenFileName(parent, paramName, currentPath);
				}
				else{
					filename = QFileDialog::getSaveFileName(parent, paramName, currentPath);
				}
				if(filename != ""){
					((TreeModel*)(index.model()))->setData(index, filename);
				}
			}
		}
		else if(paramValidator->validStringValues()->size() != 0){
			RCP<const Array<std::string> > options = paramValidator->validStringValues();
			editor = new QComboBox(parent);
			for(Array<std::string>::const_iterator itr = options->begin(); itr != options->end(); ++itr){
				static_cast<QComboBox*>(editor)->addItem(QString::fromStdString(*itr));
			}
		}
		else{
			editor = new QLineEdit(parent);
		}
	}
	else if(itemType.contains(arrayId)){
		editor = getArrayEditor(index, getArrayType(itemType), parent);
	}
	else if(itemType.contains(twoDArrayId)){
    editor = getArrayEditor(index, getArrayType(itemType), parent, true);
  }

	return editor;
}

void Delegate::setEditorData(QWidget *editor, const QModelIndex &index) const{
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);
	QVariant value = index.model()->data(index);
	if(itemType == intId){
		static_cast<QSpinBox*>(editor)->setValue(value.toInt());
	}
	else if(itemType == shortId){
		static_cast<QSpinBox*>(editor)->setValue(value.toInt());
	}
	else if(itemType == doubleId){
		static_cast<QLineEdit*>(editor)->setText(value.toString());
	}
	else if(itemType == floatId){
		static_cast<QLineEdit*>(editor)->setText(value.toString());
	}
	else if(itemType == boolId){
		static_cast<QComboBox*>(editor)->setEditText(value.toString());
	}
	else if(itemType == stringId){
		RCP<const ParameterEntryValidator> validator = ((TreeModel*)(index.model()))->getValidator(index);
		if(is_null(validator) || validator->validStringValues()->size()==0)
			static_cast<QLineEdit*>(editor)->setText(value.toString());
	 	else
			static_cast<QComboBox*>(editor)->setEditText(value.toString());
	}
  else if(itemType.contains(arrayId)){
    setArrayWidgetData(editor, getArrayType(itemType), index);
  }
  else if(itemType.contains(twoDArrayId)){
    setArrayWidgetData(editor, getArrayType(itemType), index, true);
  }
}

void Delegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const{
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);
	if(itemType == intId){
		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, spinBox->value(), Qt::EditRole);
	}
	else if(itemType == shortId){
		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, (short)spinBox->value(), Qt::EditRole);
	}
	else if(itemType == doubleId){
		QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
		model->setData(index, lineEdit->text(), Qt::EditRole);
	}
	else if(itemType == floatId){
		QLineEdit *lineEdit = static_cast<QLineEdit*>(editor);
		model->setData(index, lineEdit->text(), Qt::EditRole);
	}
	else if(itemType == boolId){
		bool value = static_cast<QComboBox*>(editor)->currentText() 
      == getBoolEditorTrue(); 
		model->setData(index, value, Qt::EditRole);
	}
	else if(itemType == stringId){
		RCP<const ParameterEntryValidator> validator = 
			static_cast<const TreeModel*>(index.model())->getValidator(index);
		QString value;
		if(is_null(validator)){
			value = static_cast<QLineEdit*>(editor)->text();
		}
		else{
			value = static_cast<QComboBox*>(editor)->currentText(); 
		}
		model->setData(index, value, Qt::EditRole);
	}
  else if(itemType.contains(arrayId)){
    QVariant value = extractValueFromArray(editor, getArrayType(itemType));
    model->setData(index, value, Qt::EditRole);
  }
  else if(itemType.contains(twoDArrayId)){
    QVariant value = extractValueFromArray(editor, getArrayType(itemType), true);
    model->setData(index, value, Qt::EditRole);
  }
}

 

void Delegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &/*index*/) const{
	editor->setGeometry(option.rect);
}

QWidget* Delegate::getArrayEditor(const QModelIndex& index, QString type, QWidget *parent, bool isTwoD) const{
  TreeModel* model = (TreeModel*)index.model();
  QString name = model->data(
    index.sibling(index.row(),0),Qt::DisplayRole).toString();
  RCP<const ParameterEntryValidator> validator = 
    model->getValidator(index);
	if(type == intId){
    if(isTwoD){
      return new Int2DArrayWidget(name, type, validator, parent);
    }
    else{
      return new IntArrayWidget(name, type, validator, parent);
    }
	}
	else if(type == shortId){
    if(isTwoD){
      return new Short2DArrayWidget(name, type, validator, parent);
    }
    else{
      return new ShortArrayWidget(name, type, validator, parent);
    }
	}
	else if(type == doubleId){
    if(isTwoD){
      return new Double2DArrayWidget(name, type, validator, parent);
    }
    else{
      return new DoubleArrayWidget(name, type, validator, parent);
    }
  }
	else if(type == floatId){
    if(isTwoD){
      return new Float2DArrayWidget(name, type, validator, parent);
    }
    else{
      return new FloatArrayWidget(name, type, validator, parent);
    }
	}
	else if(type == stringId){
    if(isTwoD){
      return new String2DArrayWidget(name, type, validator, parent);
    }
    else{
      return new StringArrayWidget(name, type, validator, parent);
    }
	}
  else{
    return 0;
  }
}

void Delegate::setArrayWidgetData(QWidget* editor, QString type, const QModelIndex& index, bool isTwoD) const{
  QVariant newData = index.model()->data(index, TreeModel::getRawDataRole());
	if(type == intId){
    isTwoD ?
    ((Int2DArrayWidget*)editor)->initData(newData.value<TwoDArray<int> >())
    :
    ((IntArrayWidget*)editor)->initData(newData.value<Array<int> >());
	}
	else if(type == shortId){
    isTwoD ?
    ((Short2DArrayWidget*)editor)->initData(newData.value<TwoDArray<short> >())
    :
    ((ShortArrayWidget*)editor)->initData(newData.value<Array<short> >());
	}
	else if(type == doubleId){
    isTwoD ?
    ((Double2DArrayWidget*)editor)->initData(newData.value<TwoDArray<double> >())
    :
    ((DoubleArrayWidget*)editor)->initData(newData.value<Array<double> >());
  }
	else if(type == floatId){
    isTwoD ?
    ((Float2DArrayWidget*)editor)->initData(newData.value<TwoDArray<float> >())
    :
    ((FloatArrayWidget*)editor)->initData(newData.value<Array<float> >());
	}
	else if(type == stringId){
    isTwoD ?
    ((String2DArrayWidget*)editor)->initData(newData.value<TwoDArray<std::string> >())
    :
    ((StringArrayWidget*)editor)->initData(newData.value<Array<std::string> >());
	}
}

QVariant Delegate::extractValueFromArray(QWidget* editor, QString type, bool isTwoD) const
{
  if(type == intId){
    return (isTwoD ?
    QVariant::fromValue(((Int2DArrayWidget*)editor)->getData())
    :
    QVariant::fromValue(((IntArrayWidget*)editor)->getData()));
  }
  else if(type == shortId){
    return (isTwoD ?
    QVariant::fromValue(((Short2DArrayWidget*)editor)->getData())
    :
    QVariant::fromValue(((ShortArrayWidget*)editor)->getData()));
  }
  else if(type == doubleId){
    return (isTwoD ?
    QVariant::fromValue(((Double2DArrayWidget*)editor)->getData())
    :
    QVariant::fromValue(((DoubleArrayWidget*)editor)->getData()));
  }
  else if(type == floatId){
    return (isTwoD ?
    QVariant::fromValue(((Float2DArrayWidget*)editor)->getData())
    :
    QVariant::fromValue(((FloatArrayWidget*)editor)->getData()));
  }
  else if(type == stringId){
    return (isTwoD ?
    QVariant::fromValue(((String2DArrayWidget*)editor)->getData())
    :
    QVariant::fromValue(((StringArrayWidget*)editor)->getData()));
  }
  else{
    return QVariant();
  }
}

} //End namespace

