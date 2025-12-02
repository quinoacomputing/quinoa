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
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

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

/** @name Functions that compute physics variables from the numerical solution for MultiMat */
///@{

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

namespace multimat {

//! Compute bulk density for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Bulk density ready to be output to file
static tk::GetVarFn::result_type
bulkDensityOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator+=;
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
  auto r = U.extract_comp( densityDofIdx(nmat,0,rdof,0) );
  for (std::size_t k=1; k<nmat; ++k)
    r += U.extract_comp( densityDofIdx(nmat,k,rdof,0) );
  return r;
}

//! Compute bulk pressure for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Bulk pressure ready to be output to file
static tk::GetVarFn::result_type
bulkPressureOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator+=;
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
  auto p = U.extract_comp( pressureDofIdx(nmat,0,rdof,0) );
  for (std::size_t k=1; k<nmat; ++k)
    p += U.extract_comp( pressureDofIdx(nmat,k,rdof,0) );
  return p;
}

//! Compute bulk specific total energy (energy per unit mass) for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Bulk specific total energy ready to be output to file
static tk::GetVarFn::result_type
bulkSpecificTotalEnergyOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator+=;
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
  auto e = U.extract_comp( energyDofIdx(nmat,0,rdof,0) );
  for (std::size_t k=1; k<nmat; ++k)
    e += U.extract_comp( energyDofIdx(nmat,k,rdof,0) );
  return e;
}

//! Compute velocity component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Velocity component ready to be output to file
template< tk::ncomp_t dir >
tk::GetVarFn::result_type
velocityOutVar( const tk::Fields& U, std::size_t rdof )
{
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
  return U.extract_comp( velocityDofIdx(nmat,dir,rdof,0) );
}

//! Compute material indicator function for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Material indicator function ready to be output to file
static tk::GetVarFn::result_type
matIndicatorOutVar( const tk::Fields& U, std::size_t rdof )
{
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
  std::vector< tk::real > m(U.nunk(), 0.0);
  for (std::size_t i=0; i<U.nunk(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      m[i] += U(i, volfracDofIdx(nmat,k,rdof,0)) *
        static_cast< tk::real >(k+1);
  }
  return m;
}

//! Compute Cauchy stress component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam idir Physical direction, encoded as 0:x, 1:y, 2:z
//! \tparam jdir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Cauchy stress component ready to be output to file
template< tk::ncomp_t idir, tk::ncomp_t jdir >
tk::GetVarFn::result_type
stressOutVar( const tk::Fields& U, std::size_t rdof )
{
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

  std::vector< tk::real > cs(U.nunk(), 0.0);
  for (std::size_t e=0; e<cs.size(); ++e) {
    for (std::size_t k=0; k<nmat; ++k) {
      tk::real asigij(0.0);

      if (solidx[k] > 0) asigij =
        U(e, stressDofIdx(nmat,solidx[k],stressCmp[idir][jdir],rdof,0));

      if (idir == jdir)
        asigij -= U(e, pressureDofIdx(nmat,k,rdof,0));

      cs[e] += asigij;
    }
  }

  return cs;
}

//! Compute inverse deformation gradient tensor component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam idir Physical direction, encoded as 0:x, 1:y, 2:z
//! \tparam jdir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Inverse deformation gradient tensor component to be output to file
template< tk::ncomp_t idir, tk::ncomp_t jdir >
tk::GetVarFn::result_type
defGradOutVar( const tk::Fields& U, std::size_t rdof )
{
  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

  std::vector< tk::real > g(U.nunk(), 0.0);
  for (std::size_t e=0; e<g.size(); ++e) {
    for (std::size_t k=0; k<nmat; ++k) {
      tk::real agij(0.0);

      if (solidx[k] > 0) agij = U(e, volfracDofIdx(nmat,k,rdof,0)) *
        U(e, deformDofIdx(nmat,solidx[k],idir,jdir,rdof,0));

      g[e] += agij;
    }
  }

  return g;
}

} // multimat::

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} //inciter::

#endif // FieldOutput_h
