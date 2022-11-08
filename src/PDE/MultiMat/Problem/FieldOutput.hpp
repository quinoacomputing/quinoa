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

#include "Fields.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Return multi-material field names to be output to file
std::vector< std::string >
MultiMatFieldNames( std::size_t nmat );

//! Return field output going to file
std::vector< std::vector< tk::real > >
MultiMatFieldOutput(
  ncomp_t system,
  std::size_t nmat,
  const std::vector< EOS >& mat_blk,
  std::size_t nunk,
  std::size_t rdof,
  const std::vector< tk::real >& vol,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const tk::Fields& U,
  const tk::Fields& P );

//! Return time history field names to be output to file
std::vector< std::string > MultiMatHistNames();

//! Return diagnostic var names to be output to file
std::vector< std::string > MultiMatDiagNames(std::size_t nmat);

} //inciter::

#endif // FieldOutput_h
