// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
  \see       FaceData.h for more info.
*/
// *****************************************************************************

#include "Reorder.hpp"
#include "DerivedData.hpp"
#include "FaceData.hpp"

using inciter::FaceData;

FaceData::FaceData(
  const std::vector< std::size_t >& inpoel,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel )
  : m_bface( bface ), m_triinpoel( triinpoel )
// *****************************************************************************
//  Constructor: compute (element-face) data for internal and domain-boundary
//  faces
//! \param[in] inpoel Mesh connectivity with local IDs
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity with local IDs
// *****************************************************************************
{
  auto esup = tk::genEsup( inpoel, 4 );
  m_esuel = tk::genEsuelTet( inpoel, esup );
  auto nbfac = tk::sumvalsize( m_bface );
  m_nipfac = tk::genNipfac( 4, nbfac, m_esuel );
  m_inpofa = tk::genInpofaTet( m_nipfac, nbfac, inpoel, m_triinpoel, m_esuel );
  m_belem =  tk::genBelemTet( nbfac, m_inpofa, esup );
  m_esuf = tk::genEsuf( 4, m_nipfac, nbfac, m_belem, m_esuel );
  Assert( m_belem.size() == nbfac,
         "Number of boundary-elements and number of boundary-faces unequal" );
}

void
FaceData::genLocalFaceId( const std::vector< std::size_t >& inpoel )
// *****************************************************************************
//  Generate local face IDs for all faces (2 per face since face has 2 esuf)
// *****************************************************************************
{
  m_faceLocalId.resize(m_esuf.size());

  for (std::size_t f=0; f<m_esuf.size()/2; ++f) {
    std::array< std::size_t, 3 > inpofa_f{{
      m_inpofa[3*f],
      m_inpofa[3*f+1],
      m_inpofa[3*f+2] }};

    // left element is always interior
    Assert( m_esuf[2*f] > -1, "Interior element detected as -1" );
    std::size_t el = static_cast< std::size_t >(m_esuf[2*f]);
    std::array< std::size_t, 4 > inpoel_l{{
      inpoel[4*el],
      inpoel[4*el+1],
      inpoel[4*el+2],
      inpoel[4*el+3] }};
    m_faceLocalId[2*f] = tk::opposite_vertex_of_tet(inpoel_l, inpofa_f);

    if (f >= Nbfac()) {
      // internal face
      Assert( m_esuf[2*f+1] > -1, "Interior element detected as -1" );
      std::size_t er = static_cast< std::size_t >(m_esuf[2*f+1]);
      std::array< std::size_t, 4 > inpoel_r{{
        inpoel[4*er],
        inpoel[4*er+1],
        inpoel[4*er+2],
        inpoel[4*er+3] }};
      m_faceLocalId[2*f+1] = tk::opposite_vertex_of_tet(inpoel_r, inpofa_f);
    }
    else {
      // boundary face
      Assert( m_esuf[2*f+1] == -1, "Outside boundary element not -1" );
      m_faceLocalId[2*f+1] = -1;  // no element to the right-side of boundary
    }
  }
}
