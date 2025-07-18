// *****************************************************************************
/*!
  \file      src/PDE/Limiter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Limiters for discontiunous Galerkin methods
  \details   This file contains functions that provide limiter function
    calculations for maintaining monotonicity near solution discontinuities
    for the DG discretization.
*/
// *****************************************************************************
#ifndef Limiter_h
#define Limiter_h

#include "Integrate/Basis.hpp"
#include "Types.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "FunctionPrototypes.hpp"
#include "EoS/EOS.hpp"
#include "Kokkos_Core.hpp"

using execution_space = Kokkos::DefaultExecutionSpace;
using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

namespace inciter {

using ncomp_t = tk::ncomp_t;

//! Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1
void
WENO_P1( const std::vector< int >& esuel,
         tk::Fields& U );

//! Superbee limiter for DGP1
void
Superbee_P1( const std::vector< int >& esuel,
             const std::vector< std::size_t >& inpoel,
             const std::vector< std::size_t >& ndofel,
             const tk::UnsMesh::Coords& coord,
             tk::Fields& U );

//! Superbee limiter for multi-material DGP1
void
SuperbeeMultiMat_P1(
  const std::vector< int >& esuel,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat );

//! Kuzmin's vertex-based limiter for transport DGP1
void
VertexBasedTransport_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U );

//! Kuzmin's vertex-based limiter for single-material DGP1
void
VertexBasedCompflow_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for single-material DGP2
void
VertexBasedCompflow_P2(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  const std::vector< std::vector<tk::real> >& mtInv,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for multi-material DGP1
void
VertexBasedMultiMat_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for multi-material DGP2
void
VertexBasedMultiMat_P2(
  const bool pref,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  const std::vector< std::vector<tk::real> >& pNodalExtrm,
  const std::vector< std::vector<tk::real> >& mtInv,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for multi-material FV
void
VertexBasedMultiMat_FV(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  std::size_t nelem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< int >& srcFlag,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat );

//! Kuzmin's vertex-based limiter for multi-species DGP1
void
VertexBasedMultiSpecies_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< inciter::EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nspec,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for multi-species DGP2
void
VertexBasedMultiSpecies_P2(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nspec,
  std::vector< std::size_t >& shockmarker );

//! WENO limiter function calculation for P1 dofs
void
WENOLimiting( const tk::Fields& U,
              const std::vector< int >& esuel,
              std::size_t e,
              inciter::ncomp_t c,
              std::size_t rdof,
              tk::real cweight,
              std::array< std::vector< tk::real >, 3 >& limU );

//! Superbee limiter function calculation for P1 dofs
std::vector< tk::real >
SuperbeeLimiting( const tk::Fields& U,
                  const std::vector< int >& esuel,
                  const std::vector< std::size_t >& inpoel,
                  const tk::UnsMesh::Coords& coord,
                  std::size_t e,
                  std::size_t ndof,
                  std::size_t rdof,
                  std::size_t dof_el,
                  inciter:: ncomp_t ncomp,
                  tk::real beta_lim );

//! Kuzmin's vertex-based limiter function calculation for P1 dofs
void
VertexBasedLimiting(
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  std::size_t e,
  std::size_t rdof,
  std::size_t ,
  std::size_t ncomp,
  std::vector< tk::real >& phi,
  const std::vector< std::size_t >& VarList );

//! Kuzmin's vertex-based limiter function calculation for P2 dofs
void
VertexBasedLimiting_P2(
  const std::vector< std::vector< tk::real > >& unk,
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  std::size_t e,
  std::size_t rdof,
  std::size_t dof_el,
  std::size_t ncomp,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& NodalExtrm,
  const std::vector< std::size_t >& VarList,
  std::vector< tk::real >& phi );

//! Consistent limiter modifications for P1 dofs
void consistentMultiMatLimiting_P1( const std::size_t nmat,
  const std::size_t rdof,
  const std::size_t e,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::vector< tk::real >& phic_p1,
  std::vector< tk::real >& phic_p2 );

//! Bound preserving limiter for the volume fractions
void BoundPreservingLimiting( std::size_t nmat,
                              std::size_t ndof,
                              std::size_t e,
                              const std::vector< std::size_t >& inpoel,
                              const tk::UnsMesh::Coords& coord,
                              const tk::Fields& U,
                              std::vector< tk::real >& phic_p1,
                              std::vector< tk::real >& phic_p2 );

//! Bound preserving limiter function for the volume fractions
tk::real
BoundPreservingLimitingFunction( const tk::real min,
                                 const tk::real max,
                                 const tk::real al_gp,
                                 const tk::real al_avg );

//! Positivity preserving limiter for multi-material and multispecies solver
void PositivityLimiting( std::size_t nmat,
                         std::size_t nspec,
                         const std::vector< EOS >& mat_blk,
                         std::size_t rdof,
                         std::size_t ndof_el,
                         const std::vector< std::size_t >& ndofel,
                         std::size_t e,
                         const std::vector< std::size_t >& inpoel,
                         const tk::UnsMesh::Coords& coord,
                         const std::vector< int >& esuel,
                         const tk::Fields& U,
                         const tk::Fields& P,
                         std::vector< tk::real >& phic_p1,
                         std::vector< tk::real >& phic_p2,
                         std::vector< tk::real >& phip_p1,
                         std::vector< tk::real >& phip_p2 );

//! Positivity bounds for multi-material PDE system
void PositivityBoundsMultiMat(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t rdof,
  std::size_t e,
  const tk::Fields& U,
  const tk::Fields& P,
  const std::vector< tk::real >& state,
  const std::vector< tk::real >& sprim,
  std::vector< tk::real >& phic_bounds,
  std::vector< tk::real >& phip_bounds );

//! Positivity bounds for multispecies PDE system
void PositivityBoundsMultiSpecies(
  std::size_t nspec,
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t rdof,
  std::size_t e,
  const tk::Fields& U,
  const tk::Fields& P,
  const std::vector< tk::real >& state,
  const std::vector< tk::real >& sprim,
  std::vector< tk::real >& phic_bound,
  std::vector< tk::real >& phip_bound );

//! Positivity preserving limiter for the FV multi-material solver
void PositivityPreservingMultiMat_FV(
  const std::vector< std::size_t >& inpoel,
  std::size_t nelem,
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoFace,
  tk::Fields& U,
  tk::Fields& P );

//! Positivity preserving limiter function
tk::real
PositivityFunction( const tk::real min,
                    const tk::real u_gp,
                    const tk::real u_avg );

//! Interface indicator function, which checks element for material interface
bool
interfaceIndicator( std::size_t nmat,
  const std::vector< tk::real >& al,
  std::vector< std::size_t >& matInt );

KOKKOS_FUNCTION
bool interfaceIndicator( std::size_t nmat,
  Kokkos::View<tk::real*, memory_space> al,
  Kokkos::View<size_t*, memory_space> matInt );

//! Mark the cells that contain discontinuity according to the interface
void MarkShockCells ( const bool pref,
                      const std::size_t nelem,
                      const std::size_t nmat,
                      const std::size_t ndof,
                      const std::size_t rdof,
                      const std::vector< EOS >& mat_blk,
                      const std::vector< std::size_t >& ndofel,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      const inciter::FaceData& fd,
                      const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const tk::FluxFn& flux,
                      const std::vector< std::size_t >& solidx,
                      const tk::Fields& U,
                      const tk::Fields& P,
                      const std::set< std::size_t >& vars,
                      std::vector< std::size_t >& shockmarker );

//! Update the conservative quantities after limiting for multi-material systems
void
correctLimConservMultiMat(
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  std::size_t nmat,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  const tk::Fields& prim,
  tk::Fields& unk );

//! Update the conservative quantities after limiting for multispecies systems
void
correctLimConservMultiSpecies(
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  std::size_t nspec,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  const tk::Fields& prim,
  tk::Fields& unk );

//! Constrain material partial pressure (alpha_k * p_k)
tk::real
constrain_pressure( const std::vector< EOS >& mat_blk,
  tk::real apr,
  tk::real arho,
  tk::real alpha,
  std::size_t imat );

} // inciter::

#endif // Limiter_h
