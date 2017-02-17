
//@HEADER
// ***********************************************************************
// 
//                     New_Package Example Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef _NEWP_JAMBO_H_
#define _NEWP_JAMBO_H_
#include "New_Package_ConfigDefs.h"
#include "Newp_Jambo.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
//! Newp_Jambo: A sample class 

/*! The Newp_Jambo class prints out a "Hello World" type message using the only word of 
Swahili that I know: Jambo, which means Hello.  The sole purpose of this code is to demonstrate
the conditional compilation and linkage features of the AUTOTOOLS package of tools.

<b>A typical heading</b>
<ul>
  <li> A typical first list entry
  <li> A typical second list entry
</ul>

<b>Another typical heading</b>

*/

//=========================================================================
class Newp_Jambo {

  public:

  //@{ \name Constructors/destructors.
  //! Basic Newp_Jambo constuctor.
  /*! Creates a Newp_Jambo object and fills with default values.  

    \warning Newp_Hello is not fully translated into Swahili.

    \param Comm In 
           An Epetra Communicator 

    \return  Newp_Jambo object

  */
  Newp_Jambo(const Epetra_Comm& Comm);

  //! Newp_Jambo copy constructor.
  
  Newp_Jambo(const Newp_Jambo& Source);
  
  //@}
  
  //@{ \name Print methods


  //! Print method
  virtual void Print(ostream & os) const;
  //@}


 private:

  const Epetra_Comm& Comm_ ; 

};

#endif /* _NEWP_JAMBO_H_ */
