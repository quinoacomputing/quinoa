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
#include "FunctionPrototypes.hpp"

namespace inciter {

//! Return a map that associates user-specified strings to functions
std::map< std::string, tk::GetVarFn > CompFlowOutVarFn();

//! Return surface field names to be output to file
std::vector< std::string > CompFlowSurfNames();

//! Return surface field output going to file
std::vector< std::vector< tk::real > >
CompFlowSurfOutput( const std::vector< EOS >& mat_blk,
                    const std::map< int, std::vector< std::size_t > >& bnd,
                    const tk::Fields& U );

//! Return element surface field output (on triangle faces) going to file
std::vector< std::vector< tk::real > >
CompFlowElemSurfOutput(
  const std::vector< EOS >& mat_blk,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const tk::Fields& U );

//! Return time history field names to be output to file
std::vector< std::string > CompFlowHistNames();

//! Return time history field output evaluated at time history points
std::vector< std::vector< tk::real > >
CompFlowHistOutput( const std::vector< EOS >& mat_blk,
                    const std::vector< HistData >& h,
                    const std::vector< std::size_t >& inpoel,
                    const tk::Fields& U );

} //inciter::

#endif // CompFlowFieldOutput_h
