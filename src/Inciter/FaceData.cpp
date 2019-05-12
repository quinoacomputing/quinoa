// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
