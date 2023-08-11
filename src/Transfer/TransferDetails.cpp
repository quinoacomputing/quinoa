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

namespace exam2m {
extern CollideHandle collideHandle;
extern CProxy_M2MTransfer m2mtransferProxy;
}

using exam2m::TransferDetails;

TransferDetails::TransferDetails( CkArrayID p, MeshData d, CkCallback cb ) :
    m_firstchunk(d.m_firstchunk)
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
  //background();

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
    if (intet(colls[i].point, colls[i].source_index, N)) {
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

bool
TransferDetails::intet(const CkVector3d &point,
      std::size_t e,
      std::array< real, 4 >& N) const
  // *****************************************************************************
  //  Determine if a point is in a tetrahedron and evaluate the shapefunction
  //! \param[in] point Point coordinates
  //! \param[in] e Mesh cell index
  //! \param[in,out] N Shapefunctions evaluated at the point
  //! \return True if ppoint is in mesh cell
  //! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
  // *****************************************************************************
{
  const std::vector< std::size_t >& inpoel = *m_inpoel;
  const tk::UnsMesh::Coords& coord = *m_coord;

  // Tetrahedron node indices
  const auto A = inpoel[e*4+0];
  const auto B = inpoel[e*4+1];
  const auto C = inpoel[e*4+2];
  const auto D = inpoel[e*4+3];

  // Tetrahedron node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Point coordinates
  const auto& xp = point.x;
  const auto& yp = point.y;
  const auto& zp = point.z;

  // Evaluate linear shapefunctions at point locations using Cramer's Rule
  //    | xp |   | x1 x2 x3 x4 |   | N1 |
  //    | yp | = | y1 y2 y3 y4 | • | N2 |
  //    | zp |   | z1 z2 z3 z4 |   | N3 |
  //    | 1  |   | 1  1  1  1  |   | N4 |

  real DetX = (y[B]*z[C] - y[C]*z[B] - y[B]*z[D] + y[D]*z[B] +
      y[C]*z[D] - y[D]*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
    x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
    x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[C]*y[A]*z[D] + x[C]*y[D]*z[A] +
    x[D]*y[A]*z[C] - x[D]*y[C]*z[A] - x[B]*y[C]*z[D] + x[B]*y[D]*z[C] +
    x[C]*y[B]*z[D] - x[C]*y[D]*z[B] - x[D]*y[B]*z[C] + x[D]*y[C]*z[B];

  real DetX1 = (y[D]*z[C] - y[C]*z[D] + y[C]*zp - yp*z[C] -
      y[D]*zp + yp*z[D])*x[B] + x[C]*y[B]*z[D] - x[C]*y[D]*z[B] -
    x[D]*y[B]*z[C] + x[D]*y[C]*z[B] - x[C]*y[B]*zp + x[C]*yp*z[B] +
    xp*y[B]*z[C] - xp*y[C]*z[B] + x[D]*y[B]*zp - x[D]*yp*z[B] -
    xp*y[B]*z[D] + xp*y[D]*z[B] + x[C]*y[D]*zp - x[C]*yp*z[D] -
    x[D]*y[C]*zp + x[D]*yp*z[C] + xp*y[C]*z[D] - xp*y[D]*z[C];

  real DetX2 = (y[C]*z[D] - y[D]*z[C] - y[C]*zp + yp*z[C] +
      y[D]*zp - yp*z[D])*x[A] + x[C]*y[D]*z[A] - x[C]*y[A]*z[D] +
    x[D]*y[A]*z[C] - x[D]*y[C]*z[A] + x[C]*y[A]*zp - x[C]*yp*z[A] -
    xp*y[A]*z[C] + xp*y[C]*z[A] - x[D]*y[A]*zp + x[D]*yp*z[A] +
    xp*y[A]*z[D] - xp*y[D]*z[A] - x[C]*y[D]*zp + x[C]*yp*z[D] +
    x[D]*y[C]*zp - x[D]*yp*z[C] - xp*y[C]*z[D] + xp*y[D]*z[C];

  real DetX3 = (y[D]*z[B] - y[B]*z[D] + y[B]*zp - yp*z[B] -
      y[D]*zp + yp*z[D])*x[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
    x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[B]*y[A]*zp + x[B]*yp*z[A] +
    xp*y[A]*z[B] - xp*y[B]*z[A] + x[D]*y[A]*zp - x[D]*yp*z[A] -
    xp*y[A]*z[D] + xp*y[D]*z[A] + x[B]*y[D]*zp - x[B]*yp*z[D] -
    x[D]*y[B]*zp + x[D]*yp*z[B] + xp*y[B]*z[D] - xp*y[D]*z[B];

  real DetX4 = (y[B]*z[C] - y[C]*z[B] - y[B]*zp + yp*z[B] +
      y[C]*zp - yp*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
    x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*zp - x[B]*yp*z[A] -
    xp*y[A]*z[B] + xp*y[B]*z[A] - x[C]*y[A]*zp + x[C]*yp*z[A] +
    xp*y[A]*z[C] - xp*y[C]*z[A] - x[B]*y[C]*zp + x[B]*yp*z[C] +
    x[C]*y[B]*zp - x[C]*yp*z[B] - xp*y[B]*z[C] + xp*y[C]*z[B];

  // Shape functions evaluated at point
  N[0] = DetX1/DetX;
  N[1] = DetX2/DetX;
  N[2] = DetX3/DetX;
  N[3] = DetX4/DetX;

  // if min( N^i, 1-N^i ) > 0 for all i, point is in cell
  if ( std::min(N[0],1.0-N[0]) > 0 && std::min(N[1],1.0-N[1]) > 0 &&
      std::min(N[2],1.0-N[2]) > 0 && std::min(N[3],1.0-N[3]) > 0 )
  {
    return true;
  } else {
    return false;
  }
}

#include "NoWarning/transferdetails.def.h"
