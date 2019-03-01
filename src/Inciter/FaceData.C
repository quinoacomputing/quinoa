// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
  \see       FaceData.h for more info.
*/
// *****************************************************************************

#include <numeric>

#include "Reorder.h"
#include "DerivedData.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "FaceData.h"
#include "Centering.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::FaceData;

FaceData::FaceData(
  const std::vector< std::size_t >& inpoel,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel )
  : m_bface( bface ), m_triinpoel( triinpoel )
// *****************************************************************************
//  Constructor: compute (element-face) data on domain boundary
//! \param[in] inpoel Mesh connectivity with local IDs
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  if (ctr::Scheme().centering(scheme) == tk::Centering::ELEM) {
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
}
