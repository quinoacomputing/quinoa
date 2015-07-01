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
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <QSize>
#include "Optika_treeitem.hpp"
#include "Optika_treemodel.hpp"

namespace Optika{

TreeItem::TreeItem(const QString& name, RCP<ParameterEntry> parameter, TreeItem *parent, bool isHeader):
  name(name),
	parentItem(parent),
	parameterEntry(parameter),
  isHeader(isHeader)
{
  myTypeId = getTypeId(parameter);
  if(isHeader){
    return;
  }
	else if(myTypeId == unrecognizedId && nonnull(parameterEntry)){
		this->docString = "Sorry, but we don't recognize the type of the " + name + " parameter.\n"
		 + "No worries though. Everything should be fine.\n"
		 "We'll just go ahead and set this parameter to its default value for you."
		 "\n\nActual Documentation:\n" + QString::fromStdString(parameter->docString());
	}
	else if(nonnull(parameter)){
		this->docString = QString::fromStdString(parameter->docString());
	}
	else{
		this->docString = "";
	}
}

TreeItem::~TreeItem(){
	qDeleteAll(childItems);
}

void TreeItem::printOut() const{
	std::cout << name.toStdString() <<  ":     ";
	for(int i=0; i<childItems.size(); ++i){
		childItems.at(i)->printOut();
	}
}

void TreeItem::appendChild(TreeItem *item){
	childItems.append(item);
}

TreeItem *TreeItem::child(int row){
	return childItems.value(row);
}

int TreeItem::childCount() const{
	return childItems.count();
}

const QList<TreeItem*> TreeItem::getChildItems(){
	return childItems;
}

int TreeItem::columnCount() const{
  return 3;
}

QVariant TreeItem::data(int column, int role) const{
	if(role == Qt::ToolTipRole){
		if(name.compare(QString("Kurtis is awesome!"), Qt::CaseInsensitive) == 0){
			return QString("I know! I think I'm awesome too!\n"
			"You're pretty awesome yourself! You should send\n"
			"me an e-mail letting me know you found the easter egg.\n"
			"I'd enjoy that.\n"
			"kob0724@gmail.com or klnusbaum@gmail.com");
		}
		else if(name.compare(QString("Jim is awesome!"), Qt::CaseInsensitive) == 0){
			return QString("I know! I think he's awesome too!\n"
			"You're pretty awesome yourself! You should send\n"
			"Jim an e-mail letting him know you think he's awesome.\n"
			"He'd enjoy that.\n"
			"Tell him Kurtis sent you. jmwille@sandia.gov");
		}
		else if(name.compare(QString("Dr. Heroux is awesome!"), Qt::CaseInsensitive) == 0){
			return QString("I know! I think he's awesome too!\n"
			"You're pretty awesome yourself! You should send\n"
			"Dr. Heroux an e-mail letting him know you think he's awesome.\n"
			"He'd enjoy that.\n"
			"Tell him Kurtis sent you. maherou@sandia.gov");
		}
		return docString;
	}
	else if(role == Qt::DisplayRole && isHeader){
		if(column == 0){
			return "Parameter";
		}
		else if (column == 1){
			return "Value";
		}
		else if(column == 2){
			return "Type";
		}
  }
	else if(role == Qt::DisplayRole && myTypeId == unrecognizedId){
		if(column == 0){
			return name;
		}
		else if (column == 1){
			return QVariant("N/A");
		}
		else if(column == 2){
			return QVariant(unrecognizedId);
		}
	}
	else if(role == Qt::DisplayRole){
    if(column == 0){
      return name;
    }
    else if(column == 1 && 
      nonnull(parameterEntry) && 
      parameterEntry->isTwoDArray()
    )
    {
      return QString("Click to view 2D Array");
    }
    else if(column == 1 && 
      nonnull(parameterEntry) &&
      !parameterEntry->isList()
    )
    {
      std::string str = toString(parameterEntry->getAny());
      if(parameterEntry->isType<bool>()) {
        if(str == "0")
          str = "false";
        else if(str == "1")
          str = "true";
      }
      return QString::fromStdString(str);
    }
    else if(column == 2){
      return myTypeId;
    }
	}
  else if(role == TreeModel::getRawDataRole()){
    if(column == 1 && nonnull(parameterEntry) && parameterEntry->isArray()){
      return arrayEntryToVariant(parameterEntry, 
        getArrayType(myTypeId));
    }
    else if(column == 1 && nonnull(parameterEntry) && parameterEntry->isTwoDArray()){
      return arrayEntryToVariant(parameterEntry,
        getArrayType(myTypeId), true);
    }
    /*else{
      return parameterEntry->getAny();
    }*/
  }
	return QVariant();
}

TreeItem* TreeItem::parent(){
	return parentItem;
}

int TreeItem::row() const{
	if(parentItem){
		return parentItem->childItems.indexOf(const_cast<TreeItem*>(this));
	}
	return 0;
}

bool TreeItem::hasValidValue() const{
	if(is_null(parameterEntry->validator())){
		return true;
  }
	else{
		try{
			parameterEntry->validator()->validate(*parameterEntry, data(0).toString().toStdString(),
							      parentItem->data(0,Qt::DisplayRole).toString().toStdString());
			return true;
		}
		catch(std::exception& /*e*/){
			return false;
		}
	}
	//should never get here
	return true;

}

QString TreeItem::getCurrentInvalidValueMessage() const{
  if(parameterEntry->validator() == null){
    return "";
  }
	try{
		parameterEntry->validator()->validate(*parameterEntry, data(0).toString().toStdString(),
						      parentItem->data(0,Qt::DisplayRole).toString().toStdString());
    return "";
	}
	catch(std::exception& e){
		return QString::fromStdString(e.what());
	}
}
	

bool TreeItem::changeValue(QVariant value){
	if(myTypeId == intId){
    int newValue = value.value<int>();
    if(newValue != getValue<int>(*parameterEntry)){
		  parameterEntry->setValue(newValue, false, parameterEntry->docString(), parameterEntry->validator());
    }
	}
	else if(myTypeId == shortId){
    short newValue = value.value<short>();
    if(newValue != getValue<short>(*parameterEntry)){
		  parameterEntry->setValue(newValue, false, parameterEntry->docString(), parameterEntry->validator());
    }
	}
	else if(myTypeId == doubleId){
    double newValue = value.value<double>();
    if(newValue != getValue<double>(*parameterEntry)){
		  parameterEntry->setValue(newValue, false, parameterEntry->docString(), parameterEntry->validator());
    }
	}
	else if(myTypeId == floatId){
    float newValue = value.value<float>();
    if(newValue != getValue<float>(*parameterEntry)){
		  parameterEntry->setValue(newValue, false, parameterEntry->docString(), parameterEntry->validator());
    }
	}
	else if(myTypeId == boolId){
    bool newValue = value.value<bool>();
    if(newValue != getValue<bool>(*parameterEntry)){
		  parameterEntry->setValue(newValue, false, parameterEntry->docString(), parameterEntry->validator());
    }
	}
	else if(myTypeId == stringId){
    std::string newValue = value.toString().toStdString();
    if(newValue != getValue<std::string>(*parameterEntry)){
		  parameterEntry->setValue(newValue, false, parameterEntry->docString(), parameterEntry->validator());
    }
	}
	else if(myTypeId.contains(arrayId)){
		changeValueForArray(value, getArrayType(myTypeId));
	}
	else if(myTypeId.contains(twoDArrayId)){
		changeValueForArray(value, getArrayType(myTypeId), true);
	}

	return true;
}

void TreeItem::setValidator(RCP<const ParameterEntryValidator> validator){
	parameterEntry->setValidator(validator);
}

void TreeItem::changeValueForArray(QVariant value, QString type, bool twoD){
	if(type == intId){
    twoD ? 
		parameterEntry->setValue(value.value<TwoDArray<int> >(), false,
					 parameterEntry->docString(), parameterEntry->validator())
    :
		parameterEntry->setValue(value.value<Array<int> >(), false,
					 parameterEntry->docString(), parameterEntry->validator());
	}
	else if(type == shortId){
    twoD ? 
		parameterEntry->setValue(value.value<TwoDArray<short> >(), false,
					 parameterEntry->docString(), parameterEntry->validator())
    :
		parameterEntry->setValue(value.value<Array<short> >(), false,
					 parameterEntry->docString(), parameterEntry->validator());
	}
	else if(type == doubleId){
    twoD ? 
		parameterEntry->setValue(value.value<TwoDArray<double> >(), false,
					 parameterEntry->docString(), parameterEntry->validator())
    :
		parameterEntry->setValue(value.value<Array<double> >(), false,
					 parameterEntry->docString(), parameterEntry->validator());
	}
	else if(type == floatId){
    twoD ? 
		parameterEntry->setValue(value.value<TwoDArray<float> >(), false,
					 parameterEntry->docString(), parameterEntry->validator())
    :
		parameterEntry->setValue(value.value<Array<float> >(), false,
					 parameterEntry->docString(), parameterEntry->validator());
	}
	else if(type == stringId){
    twoD ? 
		parameterEntry->setValue(value.value<TwoDArray<std::string> >(), false,
					 parameterEntry->docString(), parameterEntry->validator())
    :
		parameterEntry->setValue(value.value<Array<std::string> >(), false,
					 parameterEntry->docString(), parameterEntry->validator());
	}
}

QString TreeItem::getTypeId(const RCP<const ParameterEntry> parameter){
  if(parameter.is_null()){
    return unrecognizedId;
  }
  else if(parameter->isList()){
    return listId;
  }
	else if(parameter->isType<int>()){
    return intId;
	}
	else if(parameter->isType<short>()){
		return shortId;
	}
	else if(parameter->isType<double>()){
		return doubleId;
	}
	else if(parameter->isType<float>()){
		return floatId;
	}
	else if(parameter->isType<bool>()){
		return boolId;
	}
	else if(parameter->isType<std::string>()){
		return stringId;
	}
	else if(parameter->isArray()){
		QString determinedId = determineArrayType(parameter);
		if( determinedId != unrecognizedId){
			return QString(arrayId + " "+ determinedId);
		}
		else{
			return unrecognizedId;
		}
	}
  else if(parameter->isTwoDArray()){
		QString determinedId = determineArrayType(parameter, true);
		if(determinedId != unrecognizedId){
			return QString(twoDArrayId + " "+ determinedId);
		}
		else{
			return unrecognizedId;
		}
  }
	else{
	  return unrecognizedId;
	}
  //Should never get here
  //This is here to avoid compiler warnings
  return unrecognizedId;
}


} //end namespace
