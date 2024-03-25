// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for multi-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible multi-material equations.
*/
// *****************************************************************************
#ifndef FieldOutput_h
#define FieldOutput_h

#include "Types.hpp"
#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "FaceData.hpp"
#include "FunctionPrototypes.hpp"

namespace inciter {

using ncomp_t = tk::ncomp_t;

//! Return a map that associates user-specified strings to functions
std::map< std::string, tk::GetVarFn > MultiMatOutVarFn();

//! Return multi-material field names to be output to file
std::vector< std::string >
MultiMatFieldNames( std::size_t nmat );

//! Return field output going to file
std::vector< std::vector< tk::real > >
MultiMatFieldOutput(
  ncomp_t,
  std::size_t nmat,
  const std::vector< EOS >& mat_blk,
  std::size_t nunk,
  std::size_t rdof,
  const std::vector< tk::real >& vol,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const tk::Fields& U,
  const tk::Fields& P );


//! Return surface field names to be output to file
std::vector< std::string > MultiMatSurfNames();

//! Return element surface field output (on triangle faces) going to file
std::vector< std::vector< tk::real > >
MultiMatSurfOutput(
  const std::size_t nmat,
  const std::size_t rdof,
  const FaceData& fd,
  const tk::Fields& U,
  const tk::Fields& P );

//! Return time history field names to be output to file
std::vector< std::string > MultiMatHistNames();

//! Return diagnostic var names to be output to file
std::vector< std::string > MultiMatDiagNames(std::size_t nmat);

} //inciter::

#endif // FieldOutput_h
