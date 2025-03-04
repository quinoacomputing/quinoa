// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for multi-species equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible multi-species equations.
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
std::map< std::string, tk::GetVarFn > MultiSpeciesOutVarFn();

//! Return multi-species field names to be output to file
std::vector< std::string >
MultiSpeciesFieldNames( std::size_t nspec );

//! Return surface field names to be output to file
std::vector< std::string > MultiSpeciesSurfNames();

//! Return element surface field output (on triangle faces) going to file
std::vector< std::vector< tk::real > >
MultiSpeciesSurfOutput(
  const std::size_t nspec,
  const std::size_t rdof,
  const FaceData& fd,
  const tk::Fields& U,
  const tk::Fields& P );

//! Return time history field names to be output to file
std::vector< std::string > MultiSpeciesHistNames();

//! Return diagnostic var names to be output to file
std::vector< std::string > MultiSpeciesDiagNames(std::size_t nspec);

} //inciter::

#endif // FieldOutput_h
