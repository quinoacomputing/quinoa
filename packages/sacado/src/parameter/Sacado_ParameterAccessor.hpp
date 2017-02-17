// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_PARAMETERACCESSOR_HPP
#define SACADO_PARAMETERACCESSOR_HPP

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"

namespace Sacado {

  template <typename EvalType, typename EvalTypeTraits>
  class ParameterRegistration;

  /*!
   * \brief Abstract class that provides access to a parameter
   * value in a code for the parameter library. An object of this
   * type is required to construct a ParameterRegistration object.
   */
  template<typename EvalType,
           typename EvalTypeTraits = DefaultEvalTypeTraits>
  class ParameterAccessor {
  private:
    typedef typename EvalTypeTraits::template apply<EvalType>::type ScalarT;

  public:

    typedef ScalarParameterLibrary<EvalTypeTraits> ParamLib;

    virtual ~ParameterAccessor() {};

    //! Method that returns a reference to the parameter value given the name
    //! The ParameterLibrary call this method when a parameter value changes
    virtual ScalarT& getValue(const std::string &n) = 0;

    //! Method that returns a reference to the parameter value given the name
    //! The ParameterLibrary call this method when a parameter value changes
    virtual void setValue(const std::string &n, const ScalarT& v) {
      getValue(n) = v;
    }

    void registerSacadoParameter(const std::string& name,
                                 ParamLib& paramLib);

    void registerSacadoParameter(const std::string& name,
                                 const Teuchos::RCP<ParamLib>& paramLib);

  private:
    std::vector< Teuchos::RCP< ParameterRegistration<EvalType, EvalTypeTraits> > > pr_;
  };
}

// Include implementation
#include "Sacado_ParameterAccessorImp.hpp"

#endif
