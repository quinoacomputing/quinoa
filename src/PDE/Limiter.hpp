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

#include "Types.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "EoS/EoS_Base.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1
void
WENO_P1( const std::vector< int >& esuel,
         inciter::ncomp_t offset,
         tk::Fields& U );

//! Superbee limiter for DGP1
void
Superbee_P1( const std::vector< int >& esuel,
             const std::vector< std::size_t >& inpoel,
             const std::vector< std::size_t >& ndofel,
             inciter::ncomp_t offset,
             const tk::UnsMesh::Coords& coord,
             tk::Fields& U );

//! Superbee limiter for multi-material DGP1
void
SuperbeeMultiMat_P1(
  const std::vector< int >& esuel,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t system,
  inciter::ncomp_t offset,
  const tk::UnsMesh::Coords& coord,
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
  std::size_t system,
  std::size_t offset,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U );

//! Kuzmin's vertex-based limiter for single-material DGP1
void
VertexBasedCompflow_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  std::size_t offset,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for single-material DGP2
void
VertexBasedCompflow_P2(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  std::size_t offset,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  const std::vector< std::vector<tk::real> >& mtInv,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for multi-material DGP1
void
VertexBasedMultiMat_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  std::size_t offset,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat,
  std::vector< std::size_t >& shockmarker );

//! Kuzmin's vertex-based limiter for multi-material DGP2
void
VertexBasedMultiMat_P2(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  std::size_t offset,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  const std::vector< std::vector<tk::real> >& pNodalExtrm,
  const std::vector< std::vector<tk::real> >& mtInv,
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
  std::size_t system,
  std::size_t offset,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat );

//! WENO limiter function calculation for P1 dofs
void
WENOLimiting( const tk::Fields& U,
              const std::vector< int >& esuel,
              std::size_t e,
              inciter::ncomp_t c,
              std::size_t rdof,
              inciter::ncomp_t offset,
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
                  inciter::ncomp_t offset,
                  inciter:: ncomp_t ncomp,
                  tk::real beta_lim );

//! Kuzmin's vertex-based limiter function calculation for P1 dofs
void
VertexBasedLimiting( const std::vector< std::vector< tk::real > >& unk,
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  std::size_t e,
  std::size_t rdof,
  std::size_t ,
  std::size_t offset,
  std::size_t ncomp,
  std::vector< tk::real >& phi,
  const std::array< std::size_t, 2 >& VarRange );

//! Kuzmin's vertex-based limiter function calculation for P2 dofs
void
VertexBasedLimiting_P2( const std::vector< std::vector< tk::real > >& unk,
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  std::size_t e,
  std::size_t rdof,
  std::size_t dof_el,
  std::size_t offset,
  std::size_t ncomp,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& NodalExtrm,
  const std::array< std::size_t, 2 >& VarRange,
  std::vector< tk::real >& phi );

//! Consistent limiter modifications for P1 dofs
void consistentMultiMatLimiting_P1( const std::size_t nmat,
  const ncomp_t offset,
  const std::size_t rdof,
  const std::size_t e,
  tk::Fields& U,
  tk::Fields& P,
  std::vector< tk::real >& phic_p1,
  std::vector< tk::real >& phic_p2 );

//! Bound preserving limiter for the volume fractions
void BoundPreservingLimiting( std::size_t nmat,
                              ncomp_t offset,
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

//! Positivity preserving limiter for multi-material solver
void PositivityLimitingMultiMat( std::size_t nmat,
                                 std::size_t system,
                                 ncomp_t offset,
                                 std::size_t ndof,
                                 std::size_t e,
                                 const std::vector< std::size_t >& inpoel,
                                 const tk::UnsMesh::Coords& coord,
                                 const tk::Fields& U,
                                 const tk::Fields& P,
                                 std::vector< tk::real >& phic_p1,
                                 std::vector< tk::real >& phic_p2,
                                 std::vector< tk::real >& phip_p1,
                                 std::vector< tk::real >& phip_p2 );

//! Positivity preserving limiter function
tk::real
PositivityLimiting( const tk::real min,
                    const tk::real u_gp,
                    const tk::real u_avg );

//! Interface indicator function, which checks element for material interface
bool
interfaceIndicator( std::size_t nmat,
  const std::vector< tk::real >& al,
  std::vector< std::size_t >& matInt );

//! Mark the cells that contain discontinuity according to the interface
void MarkShockCells ( const std::size_t nelem,
                      const std::size_t nmat,
                      const std::size_t system,
                      const std::size_t offset,
                      const std::size_t ndof,
                      const std::size_t rdof,
                      const std::vector< inciter::EoS_Base* >& mat_blk,
                      const std::vector< std::size_t >& ndofel,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      const inciter::FaceData& fd,
                      const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const tk::Fields& U,
                      const tk::Fields& P,
                      std::vector< std::size_t >& shockmarker );

//! Function to evaluate the flux
std::vector< std::array< tk::real, 3 > >
flux( const std::size_t nmat,
      ncomp_t system,
      ncomp_t ncomp,
      const std::vector< EoS_Base* >& mat_blk,
      const std::vector< tk::real >& ugp,
      const std::vector< std::array< tk::real, 3 > >& );

//! Clean up the state of trace materials for multi-material PDE system
bool
cleanTraceMultiMat(
  std::size_t nelem,
  std::size_t system,
  const std::vector< EoS_Base* >& mat_blk,
  std::size_t offset,
  const tk::Fields& geoElem,
  std::size_t nmat,
  tk::Fields& U,
  tk::Fields& P );

//! Time step restriction for multi material cell-centered schemes
tk::real
timeStepSizeMultiMat(
  const std::vector< EoS_Base* >& mat_blk,
  const std::vector< int >& esuf,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const std::size_t nelem,
  std::size_t offset,
  std::size_t nmat,
  const tk::Fields& U,
  const tk::Fields& P );

//! Update the conservative quantities after limiting for multi-material systems
void
correctLimConservMultiMat(
  std::size_t nelem,
  const std::vector< EoS_Base* >& mat_blk,
  std::size_t nmat,
  const tk::Fields& geoElem,
  const tk::Fields& prim,
  tk::Fields& unk );

} // inciter::

#endif // Limiter_h
