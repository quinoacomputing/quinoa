/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef _IFP_LOCALPRECON_H_
#define _IFP_LOCALPRECON_H_

#include "Ifpack_config.h"

// Make the values explicit, so we can match them with FORTRAN parameters.

enum LocalPreconName
{                          // arguments
    LP_LU           =  1,  //
    LP_INVERSE      =  2,  //
    LP_SVD          =  3,  // rthresh, athresh
    LP_DIAG         = 10,  //
    LP_SOR          = 12,  // omega, iterations
    LP_SSOR         = 13,  // omega, iterations
    LP_DIAGDOM      = 15,  // 
    LP_GERSH        = 16   // alpha
};

class IFPACK_DEPRECATED ifp_LocalPrecon
{
public:
    LocalPreconName name;
    int iarg1;
    int iarg2;
    double darg1;
    double darg2;
};

#endif // _IFP_LOCALPRECON_H_
