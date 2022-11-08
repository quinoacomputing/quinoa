// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for single-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible single-material equations.
*/
// *****************************************************************************
#ifndef CompFlowFieldOutput_h
#define CompFlowFieldOutput_h

#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "History.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Return field names to be output to file
std::vector< std::string > CompFlowFieldNames();

//! Return field output going to file
std::vector< std::vector< tk::real > > 
CompFlowFieldOutput( ncomp_t system,
                     const std::vector< EOS >& mat_blk,
                     std::size_t nunk,
                     std::size_t rdof,
                     const tk::Fields& U );

//! Return surface field names to be output to file
std::vector< std::string > CompFlowSurfNames();

//! Return surface field output going to file
std::vector< std::vector< tk::real > >
CompFlowSurfOutput( ncomp_t system,
                    const std::vector< EOS >& mat_blk,
                    const std::map< int, std::vector< std::size_t > >& bnd,
                    const tk::Fields& U );

//! Return time history field names to be output to file
std::vector< std::string > CompFlowHistNames();

//! Return time history field output evaluated at time history points
std::vector< std::vector< tk::real > >
CompFlowHistOutput( ncomp_t system,
                    const std::vector< EOS >& mat_blk,
                    const std::vector< HistData >& h,
                    const std::vector< std::size_t >& inpoel,
                    const tk::Fields& U );

} //inciter::

#endif // CompFlowFieldOutput_h
