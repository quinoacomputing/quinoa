// *****************************************************************************
/*!
  \file      src/Inciter/ElemCommMap.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Communication maps for element-based discretizations
  \details   This file contains functions that set up parallel communication
    maps for element-based discretization schemes, i.e. DG and FV.
*/
// *****************************************************************************
#ifndef ElemCommMap_h
#define ElemCommMap_h

#include <vector>

#include "Types.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "DerivedData.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Local face & tet IDs associated to 3 global node IDs
//! \details This map stores tetrahedron cell faces (map key) and their
//!   associated local face ID and inner local tet id adjacent to the face
//!   (map value). A face is given by 3 global node IDs.
using FaceMap =
  std::unordered_map< tk::UnsMesh::Face,  // 3 global node IDs
                      std::array< std::size_t, 2 >, // local face & tet ID
                      tk::UnsMesh::Hash<3>,
                      tk::UnsMesh::Eq<3> >;

bool leakyAdjacency(
  const FaceData& fd,
  const tk::Fields& geoFace );

bool faceMatch(
  const FaceData& fd,
  const std::vector< size_t >& inpoel,
  const tk::UnsMesh::Coords& coord );

bool receivedChBndFaces(
  int thisIndex,
  const std::unordered_map< int, FaceMap >& bndFace,
  tk::UnsMesh::FaceSet& expChBndFace,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const tk::UnsMesh::Coords& coord );

int findchare(
  const std::unordered_map< int, FaceMap >& bndFace,
  const tk::UnsMesh::Face& t );

std::size_t nodetripletMatch(
  const FaceData& fd,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const std::vector< size_t >& inpoel,
  const std::array< std::size_t, 2 >& id,
  const tk::UnsMesh::Face& t );

void addEsuf(
  std::vector< int >& esuf,
  const std::array< std::size_t, 2 >& id,
  std::size_t ghostid );

void addEsuel(
  std::vector< int >& esuel,
  [[maybe_unused]] const std::vector< int >& esuf,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const std::vector< size_t >& inpoel,
  const std::array< std::size_t, 2 >& id,
  std::size_t ghostid,
  const tk::UnsMesh::Face& t );

void addGeoFace(
  tk::Fields& geoFace,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const tk::UnsMesh::Coords& coord,
  const tk::UnsMesh::Face& t,
  const std::array< std::size_t, 2 >& id );

} // inciter::

#endif // ElemCommMap_h
