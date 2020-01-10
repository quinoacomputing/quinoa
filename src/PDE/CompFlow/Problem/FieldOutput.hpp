// *****************************************************************************
/*!
  \file      src/PDE/Compflow/Problem/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for single-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible single-material equations.
*/
// *****************************************************************************
#ifndef CompFlowFieldOutput_h
#define CompFlowFieldOutput_h

#include "Fields.hpp"
#include "EoS/EoS.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Return field names to be output to file
std::vector< std::string > CompFlowFieldNames();

//! Return field output going to file
std::vector< std::vector< tk::real > > 
CompFlowFieldOutput( ncomp_t system,
                     ncomp_t offset,
                     tk::Fields& U );

} //inciter::

#endif // CompFlowFieldOutput_h
