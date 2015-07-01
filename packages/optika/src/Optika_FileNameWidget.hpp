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
#ifndef OPTIKA_FILENAMEWIDGET_HPP_
#define OPTIKA_FILENAMEWIDGET_HPP_
#include <QWidget>
class QLabel;

/*! \file Optika_FileNameWidget.hpp
    \brief A widget used to obtain file information
    in Optika.
*/
namespace Optika{

/**
 * \brief A small widget consisting of a button and label that allows the user
 *  to select a file through a QFileDialog. The label displays
 *  the currently selected file.
 */
class FileNameWidget : public QWidget{
	Q_OBJECT
public:
  /** \name Constructors */
  //@{

	/**
	 * \brief Constructs a FileNameWidget
	 *
	 * @param currentFileName The Filename with which the widget should be 
	 * initialized.
	 * @param parent The parent widget.
	 */
	FileNameWidget(QString currentFileName=QString(), bool mustAlreadyExist=false, QWidget* parent=0);

  //@}

  //! @name Attribute/Query Methods 
  //@{

	/**
	 * \brief Gets the current filename in the widget.
	 *
	 * @return The current filename in the widget.
	 */
	QString getCurrentFileName();

  //@}

  //! @name Setter Functions
  //@{

	/**
	 * \brief Sets the current filename in the widget.
	 *
	 * @param The name to which the widget should be set.
	 */
	void setCurrentFileName(QString newName);

  //@}

public slots:
  //! \name Public Slots
  //@{

	/**
	 * \brief Opens a QFileDialog allowing the user to select a new filename.
	 */
	void getNewFileName();

  //@}

private:
  /** \name Private Members */
  //@{
  
	/**
	 * \brief The current file name stored in the list.
	 */
	QString currentFileName;

	/**
	 * \brief The label describing the file path.
	 */
	QLabel *pathLabel;
	
	/**
	 * \brief Whether or not the file name specified must already exist.
	 */
	bool mustAlreadyExist;

  //@}
};

}

#endif //OPTIKA_FILENAMEWIDGET_HPP_
