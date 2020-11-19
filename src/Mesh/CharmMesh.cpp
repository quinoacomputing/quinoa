// *****************************************************************************
/*!
  \file      src/Mesh/CharmMesh.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Chare class declaration for holding part of a mesh
  \details   Chare class declaration for holding part of a mesh, for
             interoperating with Charm++ libraries.
*/
// *****************************************************************************

#include "CharmMesh.hpp"
#include "Reorder.hpp"
#include "DerivedData.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-parameter"
#endif

#include "Controller.hpp"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

using tk::CharmMesh;

CharmMesh::CharmMesh(
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const std::map< int, std::vector< std::size_t > >& bnode ) :
  m_inpoel( inpoel ),
  m_coord( coord ),
  m_bface( bface ),
  m_triinpoel( triinpoel ),
  m_bnode( bnode )
// *****************************************************************************
//  Constructor
//! \param[in] inpoel Vector of mesh element connectivity owned (local IDs)
//! \param[in] coord Coordinates of mesh nodes
//! \param[in] bface Boundary face lists mapped to side set ids
//! \param[in] triinpoel Triangle element connectivity
//! \param[in] bnode Boundary node lists mapped to side set ids
// *****************************************************************************
{
  Assert( !inpoel.empty(), "No elements assigned to CharmMesh chare" );
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Jacobian in input mesh to CharmMesh non-positive" );
  Assert( tk::conforming( m_inpoel, m_coord ),
          "Input mesh to CharmMesh not conforming" );
}

void
CharmMesh::transferSource()
// *****************************************************************************
//  Pass source mesh to m2m transfer library
// *****************************************************************************
{
  //exam2m::setSourceTets(thisProxy, thisIndex, &m_inpoel, &m_coord, &m_u);
}

void
CharmMesh::transferDest()
// *****************************************************************************
//  Pass Mesh Data to m2m transfer library
// *****************************************************************************
{
  //exam2m::setDestPoints(thisProxy, thisIndex, &m_coord, &m_u, CkCallback(CkIndex_MeshArray::solutionFound(), thisProxy[thisIndex]));
}

void
CharmMesh::solutionFound()
// *****************************************************************************
//  Mesh transfer complete
// *****************************************************************************
{
  //contribute( m_cbw.get< tag::solutionfound >() );
}

#include "NoWarning/charmmesh.def.h"
