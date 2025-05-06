// *****************************************************************************
/*!
  \file      src/Control/Inciter/OutVar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types and associated functions to deal with output variables
  \details   Types and associated functions to deal with output variables.
*/
// *****************************************************************************
#ifndef OutVar_h
#define OutVar_h

#include "Types.hpp"
#include "Centering.hpp"
#include "Options/PDE.hpp"

namespace inciter {
namespace ctr {

//! Output variable
struct OutVar {
  using Centering = tk::Centering;
  
  Centering centering;  //!< Centering
  std::string name;     //!< Human readable name
  std::string alias;    //!< user specified alias to the name
  tk::ncomp_t field;    //!< Field ID
  std::string varFnIdx; //!< Material-based physics label + material id

  //! Constructor: initialize all state data
  //! \param[in] n Human readable name
  //! \param[in] f Field ID
  //! \param[in] c Variable centering
  //! \param[in] vn Var function name
  explicit OutVar( Centering c = Centering::NODE,
                   const std::string& n = {},
                   const std::string& a = "",
                   tk::ncomp_t f = 0,
                   const std::string& vn = "null" ) :
    centering(c), name(n), alias(a), field(f), varFnIdx(vn) {}

  /** @name Pack/Unpack: Serialize OutVar object for Charm++ */
  ///@{
  //! Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er& p ) {
    p | centering;
    p | name;
    p | alias;
    p | field;
    p | varFnIdx;
  }
  //! Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] v OutVar object reference
  friend void operator|( PUP::er& p, OutVar& v ) { v.pup(p); }
  ///@}

  //! Query if outvar is a request for an analytic solution
  //! \return True if outvar is a request for an analytic solution
  bool analytic() const { return name.find("analytic") != std::string::npos; }

  //! Query if outvar is a request for a primitive variable
  //! \param[in] iPDE PDE type
  //! \return True if outvar should be extracted from primitive variable data
  //! \details This function returns whether the requested variable is a part
  //!   of the vector of primitive variables. This changes according to which
  //!   system of PDEs is configured.
  bool primitive( const PDEType& iPDE ) const {
    bool is_prim(false);

    if (iPDE == PDEType::MULTISPECIES) {
      if (varFnIdx.find("temperature") != std::string::npos)
      { is_prim = true; }
    }
    else {
      if (varFnIdx.find("pressure") != std::string::npos ||
          varFnIdx.find("velocity") != std::string::npos ||
          varFnIdx.find("stress") != std::string::npos )
      { is_prim = true; }
      else if ( name.length() == 2 &&
        (name.find('u') != std::string::npos ||
         name.find('U') != std::string::npos ||
         name.find('p') != std::string::npos ||
         name.find('P') != std::string::npos) )
      { is_prim = true; }
    }

    return is_prim;
  }
};

//! \brief Pack/Unpack: Namespace-scope serialize OutVar object for Charm++
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] v OutVar object reference
inline void pup( PUP::er& p, OutVar& v ) { v.pup(p); }

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

//! Operator << for writing OutVar to output streams
//! \param[in,out] os Output stream to write to
//! \param[in] outvar OutVar to write
//! \return Updated output stream
static std::ostream& operator<< ( std::ostream& os, const OutVar& outvar ) {
  if (outvar.name.empty())
    os << outvar.field+1;
  else
    os << outvar.name;
  return os;
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

} // ctr::
} // inciter::

#endif // OutVar_h
