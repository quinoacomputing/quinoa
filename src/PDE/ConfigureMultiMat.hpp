// *****************************************************************************
/*!
  \file      src/PDE/ConfigureMultiMat.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for multi-material compressible
     flow PDE
  \details   Register and compile configuration for compressible multi-material
     flow PDE.
*/
// *****************************************************************************
#ifndef ConfigureMultiMat_h
#define ConfigureMultiMat_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/PDE.hpp"
#include "PDE/MultiMat/MultiMatIndexing.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

//! Register compressible flow PDEs into PDE factory
void
registerMultiMat( DGFactory& df, std::set< ctr::PDEType >& dgt );

//! Return information on the multi-material compressible flow PDE
std::vector< std::pair< std::string, std::string > >
infoMultiMat( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt );

//! \brief Assign function that computes physics variables from the
//!   numerical solution for MultiMat
void
assignMultiMatGetVars( const std::string& name, tk::GetVarFn& f );

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
//! \param[in] offset System offset specifying the position of the MultiMat
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Bulk density ready to be output to file
static tk::GetVarFn::result_type
bulkDensityOutVar( const tk::Fields& U,
                   tk::ctr::ncomp_t offset,
                   std::size_t rdof,
                   const tk::UnsMesh::Coords&,
                   const std::vector< std::size_t >&,
                   const std::vector< tk::real >& )
{
  using tk::operator+=;
  auto sys = tk::cref_find( g_inputdeck.get< tag::sys >(), offset );
  auto nmat = g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[ sys ];
  auto r = U.extract( densityDofIdx(nmat,0,rdof,0), offset );
  for (std::size_t k=1; k<nmat; ++k)
    r += U.extract( densityDofIdx(nmat,k,rdof,0), offset );
  return r;
}

//! Compute bulk pressure for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the MultiMat
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Bulk pressure ready to be output to file
static tk::GetVarFn::result_type
bulkPressureOutVar( const tk::Fields& U,
                    tk::ctr::ncomp_t offset,
                    std::size_t rdof,
                    const tk::UnsMesh::Coords&,
                    const std::vector< std::size_t >&,
                    const std::vector< tk::real >& )
{
  using tk::operator+=;
  auto sys = tk::cref_find( g_inputdeck.get< tag::sys >(), offset );
  auto nmat = g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[ sys ];
  auto p = U.extract( pressureDofIdx(nmat,0,rdof,0), offset );
  for (std::size_t k=1; k<nmat; ++k)
    p += U.extract( pressureDofIdx(nmat,k,rdof,0), offset );
  return p;
}

//! Compute bulk specific total energy (energy per unit mass) for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the MultiMat
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Bulk specific total energy ready to be output to file
static tk::GetVarFn::result_type
bulkSpecificTotalEnergyOutVar( const tk::Fields& U,
                               tk::ctr::ncomp_t offset,
                               std::size_t rdof,
                               const tk::UnsMesh::Coords&,
                               const std::vector< std::size_t >&,
                               const std::vector< tk::real >& )
{
  using tk::operator+=;
  auto sys = tk::cref_find( g_inputdeck.get< tag::sys >(), offset );
  auto nmat = g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[ sys ];
  auto e = U.extract( energyDofIdx(nmat,0,rdof,0), offset );
  for (std::size_t k=1; k<nmat; ++k)
    e += U.extract( energyDofIdx(nmat,k,rdof,0), offset );
  return e;
}

//! Compute velocity component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the MultiMat
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Velocity component ready to be output to file
template< tk::ctr::ncomp_t dir >
tk::GetVarFn::result_type
velocityOutVar( const tk::Fields& U,
                tk::ctr::ncomp_t offset,
                std::size_t rdof,
                const tk::UnsMesh::Coords&,
                const std::vector< std::size_t >&,
                const std::vector< tk::real >& )
{
  auto sys = tk::cref_find( g_inputdeck.get< tag::sys >(), offset );
  auto nmat = g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[ sys ];
  return U.extract( velocityDofIdx(nmat,dir,rdof,0), offset );
}

//! Compute material indicator function for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] offset System offset specifying the position of the MultiMat
//!   equation system among other systems
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Material indicator function ready to be output to file
static tk::GetVarFn::result_type
matIndicatorOutVar( const tk::Fields& U,
                    tk::ctr::ncomp_t offset,
                    std::size_t rdof,
                    const tk::UnsMesh::Coords&,
                    const std::vector< std::size_t >&,
                    const std::vector< tk::real >& )
{
  auto sys = tk::cref_find( g_inputdeck.get< tag::sys >(), offset );
  auto nmat = g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[ sys ];
  std::vector< tk::real > m(U.nunk(), 0.0);
  for (std::size_t i=0; i<U.nunk(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      m[i] += U(i, volfracDofIdx(nmat,k,rdof,0), offset) *
        static_cast< tk::real >(k+1);
  }
  return m;
}

} // multimat::

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} // inciter::

#endif // ConfigureMultiMat_h
