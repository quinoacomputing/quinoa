// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef ANASAZI_EIGENSOLVER_DECL_HPP
#define ANASAZI_EIGENSOLVER_DECL_HPP

/*! \file AnasaziEigensolverDecl.hpp
    \brief Forward declaration of the virtual base class Anasazi::Eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

namespace Anasazi {

  /*! \class Eigensolver
    \brief The Eigensolver is a templated virtual base class that defines the
     basic interface that any eigensolver will support.
  
     This interface is mainly concerned with providing a set of eigensolver status method that
     can be requested from any eigensolver by an StatusTest object.
  */
  template<class ScalarType, class MV, class OP>
  class Eigensolver;
}

#endif
