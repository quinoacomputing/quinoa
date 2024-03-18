// *****************************************************************************
/*!
  \file      src/Control/Inciter/NewOutVar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types and associated functions to deal with output variables
  \details   Types and associated functions to deal with output variables.
*/
// *****************************************************************************
#ifndef NewOutVar_h
#define NewOutVar_h

#include "Centering.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"

namespace inciter {
namespace ctr {

//! Output variable
struct NewOutVar {
  
  std::size_t field;       //!< Field ID
  tk::Centering centering; //!< Centering
  std::string name;        //!< User-specified field output name
  std::size_t type;        //!< Type of OutVar: 0:evolved; 1:primitive; 2:calculated
  tk::GetVarFn getvar;     //!< Function to compute variable from numerical solution

  //! Constructor: initialize all state data
  //! \param[in] f Field ID
  //! \param[in] c Variable centering
  //! \param[in] n Human readable name
  explicit NewOutVar(
    std::size_t f = 0,
    tk::Centering c = tk::Centering::NODE,
    const std::string& n = {} ) :
  field(f), centering(c), name(n) {}

  /** @name Pack/Unpack: Serialize NewOutVar object for Charm++ */
  ///@{
  //! Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er& p ) {
    p | field;
    p | centering;
    p | name;
    getvar = assignGetVars(name);
  }
  //! Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] v NewOutVar object reference
  friend void operator|( PUP::er& p, NewOutVar& v ) { v.pup(p); }
  ///@}

  //! Query if outvar is a request for an analytic solution
  //! \return True if outvar is a request for an analytic solution
  bool analytic() const { return name.find("analytic") != std::string::npos; }

  //! Query if outvar is a request for a (multimat) primitive variable
  //! \return True if outvar should be extracted from primitive variable data
  bool primitive() const {
    bool is_prim(false);
    if (type == 1) is_prim = true;
    return is_prim;
  }

  //! Assign getvarfn
  void assignGetVar() {
    getvar = assignGetVars(name);
  }
};

//! \brief Pack/Unpack: Namespace-scope serialize NewOutVar object for Charm++
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] v NewOutVar object reference
inline void pup( PUP::er& p, NewOutVar& v ) { v.pup(p); }

} // ctr::
} // inciter::

#endif // NewOutVar_h
