// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
  \see       FaceData.h for more info.
*/
// *****************************************************************************

#include <algorithm>

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
  const std::vector< std::size_t >& ginpoel,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  std::vector< std::size_t >& triinpoel )
  : m_bface( bface ), m_bnode( bnode ), m_triinpoel( triinpoel )
// *****************************************************************************
//  Constructor
//! \param[in] ginpoel Mesh element connectivity owned (global IDs) mesh chunk
//!   this chare operates on
//! \param[in] bface Map of boundary-face lists mapped to corresponding 
//!   side set ids for this mesh chunk
//! \param[in] bnode Map of boundary-node lists mapped to corresponding 
//!   side set ids for this mesh chunk
//! \param[in] triinpoel Interconnectivity of points and boundary-face in this
//!   mesh chunk
//! \details This class is created per chare-worker in 
//!   Partitioner::createWorkers(). This is done so that Discretization will 
//!   not hold all this data unnecessarily, viz. for DiagCG, where 
//!   it's not needed. It will be computed only when DG discretization is 
//!   chosen.
// *****************************************************************************
{
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  if (ctr::Scheme().centering(scheme) == tk::Centering::ELEM) {
    auto el = tk::global2local( ginpoel );   // fills inpoel, m_gid, m_lid
    const auto& inpoel = std::get< 0 >( el );
    const auto& lid = std::get< 2 >( el );

    auto esup = tk::genEsup(inpoel,4);
    m_esuel = tk::genEsuelTet( inpoel, esup );

    // Map face connectivity from global to local ids
    for (auto& i : m_triinpoel) i = tk::cref_find(lid,i);

    auto nbfac = tk::sumvalsize( m_bface );

    m_nipfac = tk::genNipfac( 4, nbfac, m_esuel );
    m_inpofa = tk::genInpofaTet( m_nipfac, nbfac, inpoel, m_triinpoel, m_esuel );
    m_belem =  tk::genBelemTet( nbfac, m_inpofa, esup );
    m_esuf = tk::genEsuf( 4, m_nipfac, nbfac, m_belem, m_esuel );

    Assert( m_belem.size() == nbfac,
           "Number of boundary-elements and number of boundary-faces unequal" );
  }
}
