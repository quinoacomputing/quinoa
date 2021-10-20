// *****************************************************************************
/*!
  \file      src/Inciter/ElemCommMap.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Communication maps for element-based discretizations
  \details   This file contains functions that set up parallel communication
    maps for element-based discretization schemes, i.e. DG and FV.
*/
// *****************************************************************************

#include "ElemCommMap.hpp"

#include "Tags.hpp"
#include "Vector.hpp"

namespace inciter {

bool leakyAdjacency(
  const FaceData& fd,
  const tk::Fields& geoFace )
// *****************************************************************************
// Perform leak-test on chare boundary faces
//! \details This function computes a surface integral over the boundary of the
//!   faces after the face adjacency communication map is completed. A non-zero
//!   vector result indicates a leak, e.g., a hole in the partition (covered by
//!   the faces of the face adjacency communication map), which indicates an
//!   error upstream in the code that sets up the face communication data
//!   structures.
//! \note Compared to tk::leakyPartition() this function performs the leak-test
//!   on the face geometry data structure enlarged by ghost faces on this
//!   partition by computing a discrete surface integral considering the
//!   physical and chare boundary faces, which should be equal to zero for a
//!   closed domain.
//! \return True if our chare face adjacency leaks.
// *****************************************************************************
{
  // Storage for surface integral over our chunk of the adjacency
  std::array< tk::real, 3 > s{{ 0.0, 0.0, 0.0 }};

  // physical boundary faces
  for (std::size_t f=0; f<fd.Nbfac(); ++f) {
    s[0] += geoFace(f,0,0) * geoFace(f,1,0);
    s[1] += geoFace(f,0,0) * geoFace(f,2,0);
    s[2] += geoFace(f,0,0) * geoFace(f,3,0);
  }

  // chare-boundary faces
  for (std::size_t f=fd.Nipfac(); f<fd.Esuf().size()/2; ++f) {
    s[0] += geoFace(f,0,0) * geoFace(f,1,0);
    s[1] += geoFace(f,0,0) * geoFace(f,2,0);
    s[2] += geoFace(f,0,0) * geoFace(f,3,0);
  }

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;
  return std::abs(s[0]) > eps || std::abs(s[1]) > eps || std::abs(s[2]) > eps;
}

bool faceMatch(
  const FaceData& fd,
  const std::vector< size_t >& inpoel,
  const tk::UnsMesh::Coords& coord )
// *****************************************************************************
// Check if esuf of chare-boundary faces matches
//! \details This function checks each chare-boundary esuf entry for the left
//!   and right elements. Then, it tries to match all vertices of these
//!   elements. Exactly three of these vertices must match if the esuf entry
//!   has been updated correctly at chare-boundaries.
//! \return True if chare-boundary faces match.
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();
  bool match(true);

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;

  for (auto f=fd.Nipfac(); f<esuf.size()/2; ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    std::size_t count = 0;

    for (std::size_t i=0; i<4; ++i)
    {
      auto ip = inpoel[4*el+i];
      for (std::size_t j=0; j<4; ++j)
      {
        auto jp = inpoel[4*er+j];
        auto xdiff = std::abs( coord[0][ip] - coord[0][jp] );
        auto ydiff = std::abs( coord[1][ip] - coord[1][jp] );
        auto zdiff = std::abs( coord[2][ip] - coord[2][jp] );

        if ( xdiff<=eps && ydiff<=eps && zdiff<=eps ) ++count;
      }
    }

    match = (match && count == 3);
  }

  return match;
}

bool receivedChBndFaces(
  int thisIndex,
  const std::unordered_map< int, FaceMap >& bndFace,
  tk::UnsMesh::FaceSet& expChBndFace,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const tk::UnsMesh::Coords& coord )
// *****************************************************************************
// Verify that all chare-boundary faces have been received
//! \return True if all chare-boundary faces have been received
// *****************************************************************************
{
  tk::UnsMesh::FaceSet recvBndFace;

  // Collect chare-boundary faces that have been received and expected
  for (const auto& c : bndFace)
    for (const auto& f : c.second)
      if (expChBndFace.find(f.first) != end(expChBndFace))
        recvBndFace.insert(f.first);

   // Collect info on expected but not received faces
   std::stringstream msg;
   for (const auto& f : expChBndFace)
     if (recvBndFace.find(f) == end(recvBndFace)) {
       const auto& x = coord[0];
       const auto& y = coord[1];
       const auto& z = coord[2];
       auto A = tk::cref_find( lid, f[0] );
       auto B = tk::cref_find( lid, f[1] );
       auto C = tk::cref_find( lid, f[2] );
       msg << '{' << A << ',' << B << ',' << C << "}:("
           << x[A] << ',' << y[A] << ',' << z[A] << ' '
           << x[B] << ',' << y[B] << ',' << z[B] << ' '
           << x[C] << ',' << y[C] << ',' << z[C] << ") ";
     }

  tk::destroy( expChBndFace );

  // Error out with info on missing faces
  auto s = msg.str();
  if (!s.empty()) {
    Throw( "FV chare " + std::to_string(thisIndex) +
           " missing face(s) {local node ids} (node coords): " + s );
  } else {
    return true;
  }
}

int findchare(
  const std::unordered_map< int, FaceMap >& bndFace,
  const tk::UnsMesh::Face& t )
// *****************************************************************************
// Find any chare for face (given by 3 global node IDs)
//! \param[in] t Face given by three global node IDs
//! \return Chare ID if found among any of the chares we communicate along
//!   faces with, -1 if the face cannot be found.
// *****************************************************************************
{
  for (const auto& cf : bndFace)
    // cppcheck-suppress useStlAlgorithm
    if (cf.second.find(t) != end(cf.second))
      return cf.first;
  return -1;
}

std::size_t nodetripletMatch(
  const FaceData& fd,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const std::vector< size_t >& inpoel,
  const std::array< std::size_t, 2 >& id,
  const tk::UnsMesh::Face& t )
// *****************************************************************************
// Check if entries in inpoel, inpofa and node-triplet are consistent
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] t node-triplet associated with the chare boundary face
//! \return number of nodes in inpoel that matched with t and inpofa
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  std::size_t counter = 0;
  for (std::size_t k=0; k<4; ++k)
  {
    auto el = esuf[ 2*id[0] ];
    auto ip = inpoel[ 4*static_cast< std::size_t >( el )+k ];
    Assert( el == static_cast< int >( id[1] ), "Mismatch in id and esuf" );
    for (std::size_t j=0; j<3; ++j)
    {
      auto jp = tk::cref_find( lid, t[j] );
      auto fp = inpofa[ 3*id[0]+(2-j) ];
      if (ip == jp && ip == fp) ++counter;
    }
  }

  return counter;
}

void addEsuf(
  std::vector< int >& esuf,
  const std::array< std::size_t, 2 >& id,
  std::size_t ghostid )
// *****************************************************************************
// Fill elements surrounding a face along chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] ghostid Local ID for ghost tet
//! \details This function extends and fills in the elements surrounding faces
//!   data structure (esuf) so that the left and right element id is filled
//!   in correctly on chare boundaries to contain the correct inner tet id and
//!   the local tet id for the outer (ghost) tet, both adjacent to the given
//!   chare-face boundary. Prior to this function, this data structure does not
//!   have yet face-element connectivity adjacent to chare-boundary faces, only
//!   for physical boundaries and internal faces that are not on the chare
//!   boundary (this latter purely as a result of mesh partitioning). The remote
//!   element id of the ghost is stored in a location that is local to our own
//!   esuf. The face numbering is such that esuf stores the element-face
//!   connectivity first for the physical-boundary faces, followed by that of
//!   the internal faces, followed by the chare-boundary faces. As a result,
//!   esuf can be used by physics algorithms in exactly the same way as would be
//!   used in serial. In serial, of course, this data structure is not extended
//!   at the end by the chare-boundaries.
// *****************************************************************************
{
  Assert( 2*id[0]+1 < esuf.size(), "Indexing out of esuf" );

  // put in inner tet id
  Assert( esuf[ 2*id[0] ] == -2 && esuf[ 2*id[0]+1 ] == -2, "Updating esuf at "
          "wrong location instead of chare-boundary" );
  esuf[ 2*id[0]+0 ] = static_cast< int >( id[1] );
  // put in local id for outer/ghost tet
  esuf[ 2*id[0]+1 ] = static_cast< int >( ghostid );
}

void addEsuel(
  std::vector< int >& esuel,
  [[maybe_unused]] const std::vector< int >& esuf,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const std::vector< size_t >& inpoel,
  const std::array< std::size_t, 2 >& id,
  std::size_t ghostid,
  const tk::UnsMesh::Face& t )
// *****************************************************************************
// Fill elements surrounding a element along chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] ghostid Local ID for ghost tet
//! \param[in] t node-triplet associated with the chare boundary face
//! \details This function updates the elements surrounding element (esuel) data
//    structure for the (inner) tets adjacent to the chare-boundaries. It fills
//    esuel of this inner tet with the local tet-id that has been assigned to
//    the outer ghost tet in FV::comGhost in place of the -1 before.
// *****************************************************************************
{
  std::array< tk::UnsMesh::Face, 4 > face;
  for (std::size_t f = 0; f<4; ++f)
    for (std::size_t i = 0; i<3; ++i)
      face[f][i] = inpoel[ id[1]*4 + tk::lpofa[f][i] ];

  tk::UnsMesh::Face tl{{ tk::cref_find( lid, t[0] ),
                         tk::cref_find( lid, t[1] ),
                         tk::cref_find( lid, t[2] ) }};

  std::size_t i(0), nmatch(0);
  for (const auto& f : face) {
    if (tk::UnsMesh::Eq< 3 >()( tl, f )) {
      Assert( esuel[ id[1]*4 + i ] == -1, "Incorrect boundary element found in "
             "esuel");
      esuel[ id[1]*4 + i ] = static_cast<int>(ghostid);
      ++nmatch;
      Assert( esuel[ id[1]*4 + i ] == esuf[ 2*id[0]+1 ], "Incorrect boundary "
             "element entered in esuel" );
      Assert( static_cast<int>(id[1]) == esuf[ 2*id[0]+0 ], "Boundary "
             "element entered in incorrect esuel location" );
    }
    ++i;
  }

  // ensure that exactly one face matched
  Assert( nmatch == 1, "Incorrect number of node-triplets (faces) matched for "
         "updating esuel; matching faces = "+ std::to_string(nmatch) );
}

void addGeoFace(
  tk::Fields& geoFace,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const tk::UnsMesh::Coords& coord,
  const tk::UnsMesh::Face& t,
  const std::array< std::size_t, 2 >& id )
// *****************************************************************************
// Fill face-geometry data along chare boundary
//! \param[in] t Face (given by 3 global node IDs) on the chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to face t
//! \details This function fills in the face geometry data along a chare
//!    boundary.
// *****************************************************************************
{
  // get global node IDs reversing order to get outward-pointing normal
  auto A = tk::cref_find( lid, t[2] );
  auto B = tk::cref_find( lid, t[1] );
  auto C = tk::cref_find( lid, t[0] );
  auto geochf = tk::geoFaceTri( {{coord[0][A], coord[0][B], coord[0][C]}},
                                {{coord[1][A], coord[1][B], coord[1][C]}},
                                {{coord[2][A], coord[2][B], coord[2][C]}} );

  for (std::size_t i=0; i<7; ++i)
    geoFace(id[0],i,0) = geochf(0,i,0);
}

}
// inciter::
