//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_TYPES_HPP
#define Rythmos_TYPES_HPP


#include "Rythmos_ConfigDefs.h"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_AbstractFactory.hpp"
#include "Teuchos_TimeMonitor.hpp"




namespace Rythmos {


// Using declarations from Teuchos

// 1/15/09 tscoffe:  Don't put using declarations for functions here!  Types are okay.

/** \brief . */
using Teuchos::Ptr;
/** \brief . */
using Teuchos::RCP;
/** \brief . */
using Teuchos::FancyOStream;
/** \brief . */
using Teuchos::ArrayView;
/** \brief . */
using Teuchos::Array;
/** \brief . */
using Teuchos::ArrayRCP;
/** \brief . */
using Teuchos::ParameterList;
/** \brief . */
using Teuchos::ScalarTraits;
/** \brief . */
using Teuchos::TypeNameTraits;
/** \brief . */
using Teuchos::AbstractFactory;



// Using declarations from Teuchos


/** \brief . */
using Thyra::VectorBase;


namespace Exceptions {


/** brief Base for all Rythmos exceptions that are not just simple usage
 * errors (i.e. failing preconditions).
 */
class ExceptionBase : public std::runtime_error
{public: ExceptionBase(const std::string& what_arg) : std::runtime_error(what_arg) {}};


} // namespace Exceptions


} // namespace Rythmos


#endif // Rythmos_TYPES_HPP



