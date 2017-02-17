// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
// @HEADER
//

#ifndef ANASAZI_STATUS_TEST_DECL_HPP
#define ANASAZI_STATUS_TEST_DECL_HPP

/*!
  \file AnasaziStatusTestDecl.hpp
  \brief Forward declaration of pure virtual base class Anasazi::StatusTest.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

namespace Anasazi {

  /// \class StatusTest
  /// \brief Common interface of stopping criteria for Anasazi's solvers.
  ///
  /// StatusTest is an interface that can be implemented to create
  /// convergence tests for all Anasazi solvers. Almost any kind of test
  /// can be expressed using this mechanism, including composite tests
  /// (see StatusTestCombo).
  template <class ScalarType, class MV, class OP>
  class StatusTest;
}


#endif
