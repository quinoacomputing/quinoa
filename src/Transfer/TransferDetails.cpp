// *****************************************************************************
/*!
  \file      src/Transfer/TransferDetails.cpp
  \copyright 2020 Charmworks, Inc.
             All rights reserved. See the LICENSE file for details.
  \brief     Chare class declaration for mesh transfer workers holding part of a
    mesh
  \details   Chare class declaration for mesh transfer workers holding part of a
    mesh.
*/
// *****************************************************************************

#include "TransferDetails.hpp"
#include "Reorder.hpp"
#include "DerivedData.hpp"
#include "M2MTransfer.hpp"

#include "collidecharm.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

PUPbytes(Collision);

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

extern exam2m::CProxy_M2MTransfer m2mtransferProxy;

using exam2m::TransferDetails;

TransferDetails::TransferDetails(
  CkArrayID p,
  MeshData d,
  CollideHandle ch,
  CkCallback cb ) :
    m_firstchunk(d.m_firstchunk),
    collideHandle(ch)
// *****************************************************************************
//  Constructor
//! \param[in] firstchunk Chunk ID used for the collision detection library
//! \param[in] cb Callback to inform application that the library is ready
// *****************************************************************************
{
  CollideRegister(collideHandle, m_firstchunk + thisIndex);
  d.m_proxy = thisProxy;
  m2mtransferProxy.ckLocalBranch()->setMesh( p, d );
  contribute(cb);
}

void
TransferDetails::setSourceTets(
    std::vector< std::size_t>* inpoel,
    tk::UnsMesh::Coords* coords,
    const tk::Fields& u )
// *****************************************************************************
//  Set the data for the source tetrahedrons to be collided
//! \param[in] inpoel Pointer to the connectivity data for the source mesh
//! \param[in] coords Pointer to the coordinate data for the source mesh
//! \param[in] u Pointer to the solution data for the source mesh
// *****************************************************************************
{
  m_coord = coords;
  m_u = const_cast< tk::Fields* >( &u );
  m_inpoel = inpoel;

  // Send tetrahedron data to the collision detection library
  collideTets();
}

void
TransferDetails::setDestPoints(
    tk::UnsMesh::Coords* coords,
    tk::Fields& u,
    CkCallback cb )
// *****************************************************************************
//  Set the data for the destination points to be collided
//! \param[in] coords Pointer to the coordinate data for the destination mesh
//! \param[in,out] u Pointer to the solution data for the destination mesh
//! \param[in] cb Callback to call once this chare received all solution data
// *****************************************************************************
{
  m_coord = coords;
  m_u = static_cast< tk::Fields* >( &u );
  m_donecb = cb;

  // Initialize msg counters, callback, and background solution data
  m_numsent = 0;
  m_numreceived = 0;
  background();

  // Send vertex data to the collision detection library
  collideVertices();
}

void
TransferDetails::background()
// *****************************************************************************
// Initialize dest mesh solution with background data
//! \details This is useful to see what points did not receive solution.
// *****************************************************************************
{
  tk::Fields& u = *m_u;
  for (std::size_t i = 0; i < u.nunk(); ++i) u(i,0) = -1.0;
}

void
TransferDetails::collideVertices()
// *****************************************************************************
// Pass vertex information to the collision detection library
// *****************************************************************************
{
  const tk::UnsMesh::Coords& coord = *m_coord;
  auto nVertices = coord[0].size();
  std::size_t nBoxes = 0;
  std::vector< bbox3d > boxes( nVertices );
  std::vector< int > prio( nVertices );
  auto firstchunk = static_cast< int >( m_firstchunk );
  for (std::size_t i=0; i<nVertices; ++i) {
    boxes[nBoxes].empty();
    boxes[nBoxes].add(CkVector3d(coord[0][i], coord[1][i], coord[2][i]));
    prio[nBoxes] = firstchunk;
    ++nBoxes;
  }
  CollideBoxesPrio( collideHandle, firstchunk + thisIndex,
                    static_cast<int>(nBoxes), boxes.data(), prio.data() );
}

void
TransferDetails::collideTets() const
// *****************************************************************************
// Pass tet information to the collision detection library
// *****************************************************************************
{
  const std::vector< std::size_t >& inpoel = *m_inpoel;
  const tk::UnsMesh::Coords& coord = *m_coord;
  auto nBoxes = inpoel.size() / 4;
  std::vector< bbox3d > boxes( nBoxes );
  std::vector< int > prio( nBoxes );
  auto firstchunk = static_cast< int >( m_firstchunk );
  for (std::size_t i=0; i<nBoxes; ++i) {
    boxes[i].empty();
    prio[i] = firstchunk;
    for (std::size_t j=0; j<4; ++j) {
      // Get index of the jth point of the ith tet
      auto p = inpoel[i * 4 + j];
      // Add that point to the tets bounding box
      boxes[i].add(CkVector3d(coord[0][p], coord[1][p], coord[2][p]));
    }
  }
  CollideBoxesPrio( collideHandle, firstchunk + thisIndex,
                    static_cast<int>(nBoxes), boxes.data(), prio.data() );
}

void
TransferDetails::processCollisions(
    CProxy_TransferDetails proxy,
    int numchares,
    int chunkoffset,
    int nColl,
    Collision* colls )
// *****************************************************************************
//  Process potential collisions by sending my points to the source mesh chares
//  that they potentially collide with.
//! \param[in] proxy Proxy for the source mesh chares
//! \param[in] numchares Number of chares in the source mesh chare array
//! \param[in] chunkoffset First chunk ID of the source mesh
//! \param[in] nColl Number of potential collisions to process
//! \param[in] colls List of potential collisions
// *****************************************************************************
{
  const tk::UnsMesh::Coords& coord = *m_coord;
  int mychunk = thisIndex + m_firstchunk;

  std::vector< std::vector< PotentialCollision > >
    pColls( static_cast<std::size_t>(numchares) );

  // Separate potential collisions into lists based on the source mesh chare
  // that is involved in the potential collision
  for (int i = 0; i < nColl; i++) {
    int chareindex;
    PotentialCollision pColl;
    if (colls[i].A.chunk == mychunk) {
      chareindex = colls[i].B.chunk - chunkoffset;
      pColl.dest_index = static_cast<std::size_t>(colls[i].A.number);
      pColl.source_index = static_cast<std::size_t>(colls[i].B.number);
    } else {
      chareindex = colls[i].A.chunk - chunkoffset;
      pColl.dest_index = static_cast<std::size_t>(colls[i].B.number);
      pColl.source_index = static_cast<std::size_t>(colls[i].A.number);
    }

    #if defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wdeprecated-copy"
    #endif
    pColl.point = { coord[0][pColl.dest_index],
                    coord[1][pColl.dest_index],
                    coord[2][pColl.dest_index] };
    #if defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #endif

    pColls[ static_cast<std::size_t>(chareindex) ].push_back( pColl );
  }

  // Send out the lists of potential collisions to the source mesh chares
  for (int i = 0; i < numchares; i++) {
    auto I = static_cast< std::size_t >( i );
    m_numsent++;
    proxy[i].determineActualCollisions( thisProxy,
                                        thisIndex,
                                        static_cast<int>(pColls[I].size()),
                                        pColls[I].data() );
  }
}

void
TransferDetails::determineActualCollisions(
    CProxy_TransferDetails proxy,
    int index,
    int nColls,
    PotentialCollision* colls ) const
// *****************************************************************************
//  Identify actual collisions by calling intet on all possible collisions, and
//  interpolate solution values to send back to the destination mesh.
//! \param[in] proxy The proxy of the destination mesh chare array
//! \param[in] index The index in proxy to return the solution data to
//! \param[in] nColls Number of collisions to be checked
//! \param[in] colls List of potential collisions
// *****************************************************************************
{
  const std::vector< std::size_t >& inpoel = *m_inpoel;
  tk::Fields& u = *m_u;
  //CkPrintf("Source chare %i received data for %i potential collisions\n",
  //    thisIndex, nColls);

  std::array< real, 4 > N;
  int numInTet = 0;
  std::vector< SolutionData > return_data;

  // Iterate over my potential collisions and call intet() to determine
  // if an actual collision occurred, and if so interpolate solution to dest
  for (int i = 0; i < nColls; i++) {
    std::vector< tk::real > point
      {colls[i].point.x, colls[i].point.y, colls[i].point.z};
    if (tk::intet(*m_coord, *m_inpoel, point, colls[i].source_index, N))
    {
      numInTet++;
      SolutionData data;
      data.dest_index = colls[i].dest_index;
      auto e = colls[i].source_index;
      const auto A = inpoel[e*4+0];
      const auto B = inpoel[e*4+1];
      const auto C = inpoel[e*4+2];
      const auto D = inpoel[e*4+3];
      data.solution.resize( u.nprop() );
      for (std::size_t c=0; c<u.nprop(); ++c) {
        data.solution[c] = N[0]*u(A,c) + N[1]*u(B,c) + N[2]*u(C,c) + N[3]*u(D,c);
      }
      return_data.push_back(data);
    }
  }
  //CkPrintf("Source chare %i found %i/%i actual collisions\n",
  //    thisIndex, numInTet, nColls);
  // Send the solution data for the actual collisions back to the dest mesh
  proxy[index].transferSolution( return_data );
}

void
TransferDetails::transferSolution( const std::vector< SolutionData >& soln )
// *****************************************************************************
//  Receive the solution data for destination mesh points that collided with the
//  source mesh tetrahedrons
//! \param[in] soln List of solutions
// *****************************************************************************
{
  tk::Fields& u = *m_u;

  for (std::size_t i=0; i<soln.size(); ++i) {
    for (std::size_t c=0; c<u.nprop(); ++c) {
      u(soln[i].dest_index,c) = soln[i].solution[c];
    }
  }

  // Inform the caller if we've received all solution data
  m_numreceived++;
  if (m_numreceived == m_numsent) {
    m_donecb.send();
  }
}

#include "NoWarning/transferdetails.def.h"
