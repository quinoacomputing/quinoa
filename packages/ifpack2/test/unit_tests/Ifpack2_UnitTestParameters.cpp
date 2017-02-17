// ***********************************************************************
// 
//      Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


/*! \file Ifpack2_UnitTestParameters.cpp

\brief Ifpack2 Parameters testing.

*/

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <iostream>
#include <Ifpack2_Parameters.hpp>

TEUCHOS_UNIT_TEST(Ifpack2Parameters, Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  Teuchos::ParameterList params;

  params.set("fact: iluk level-of-fill", (int) 2);

  Teuchos::ParameterList validparams;

  TEST_NOTHROW(Ifpack2::getValidParameters(validparams));

  params.validateParameters(validparams);

  int level_of_fill = 0;

  //call getParameter with a wrong name:
  Ifpack2::getParameter(params, "level-of-fill", level_of_fill);
  TEST_EQUALITY(level_of_fill, 0)

  //call getParameter with a valid name:
  Ifpack2::getParameter(params, "fact: iluk level-of-fill", level_of_fill);
  TEST_EQUALITY(level_of_fill, 2)
}

