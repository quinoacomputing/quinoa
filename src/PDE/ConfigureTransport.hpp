// *****************************************************************************
/*!
  \file      src/PDE/ConfigureTransport.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the Transport PDE
  \details   Register and compile configuration on the Transport PDE.
*/
// *****************************************************************************
#ifndef ConfigureTransport_h
#define ConfigureTransport_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/PDE.hpp"

namespace inciter {

//! Register transport PDEs into PDE factory
void
registerTransport( CGFactory& cf,
                   DGFactory& df,
                   std::set< ctr::PDEType >& cgt,
                   std::set< ctr::PDEType >& dgt );

//! Return information on the transport PDE
std::vector< std::pair< std::string, std::string > >
infoTransport( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt );

//! \brief Assign function that computes physics variables from the
//!   numerical solution for MultiMat
void
assignTransportGetVars( const std::string& name, tk::GetVarFn& f );

/** @name Functions that compute physics variables from the numerical solution for Transport */
///@{

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

namespace transport {

//! Compute material indicator function for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the Transport
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Material indicator function ready to be output to file
static tk::GetVarFn::result_type
matIndicatorOutVar( const tk::Fields& U, tk::ctr::ncomp_t offset,
                    std::size_t rdof )
{
  auto ncomp = U.nprop()/rdof;
  std::vector< tk::real > m(U.nunk(), 0.0);
  for (std::size_t i=0; i<U.nunk(); ++i) {
    for (std::size_t k=0; k<ncomp; ++k)
      m[i] += U(i, rdof*k) * static_cast< tk::real >(k+1);
  }
  return m;
}

} // transport::

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} // inciter::

#endif // ConfigureTransport_h
