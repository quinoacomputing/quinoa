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

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::FaceData;

FaceData::FaceData(
  const std::vector< std::size_t >& conn,
  const std::unordered_map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel )
  : m_bface( bface ), m_triinpoel( triinpoel )
// *****************************************************************************
//  Constructor
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//!   mesh chunk we operate on
//! \param[in] bface Map of boundary-face lists mapped to corresponding 
//!   side set ids for this mesh chunk
//! \param[in] triinpoel Interconnectivity of points and boundary-face in this
//!   mesh chunk
//! \details This class is created per chare-worker in 
//!   Partitioner::createWorkers(). This is done so that Discretization will 
//!   not hold all this data unnecessarily, viz. for MatCG and DiagCG, where 
//!   it's not needed. It will be computed only when DG discretization is 
//!   chosen.
// *****************************************************************************
{
  if (g_inputdeck.get< tag::selected, tag::scheme >() == ctr::SchemeType::DG) {

    auto el = tk::global2local( conn );   // fills inpoel, m_gid, m_lid
    auto inpoel = std::get< 0 >( el );
    auto lid = std::get< 2 >( el );

    auto esup = tk::genEsup(inpoel,4);
    m_esuel = tk::genEsuelTet( inpoel, esup );

    // Mapping m_triinpoel from global renumbered ids to local renumbered ids
    for (auto& i : m_triinpoel) i = tk::cref_find(lid,i);

    auto nbfac = numBndFaces();

    m_ntfac = tk::genNtfac( 4, nbfac, m_esuel );
    m_inpofa = tk::genInpofaTet( m_ntfac, nbfac, inpoel, m_triinpoel, m_esuel );
    m_belem =  tk::genBelemTet( nbfac, m_inpofa, esup );
    m_esuf = tk::genEsuf( 4, m_ntfac, nbfac, m_belem, m_esuel );

    Assert( m_belem.size() == nbfac,
            "Number of boundary-elements and number of boundary-faces unequal" );
  }
}

std::size_t
FaceData::numBndFaces() const
// *****************************************************************************
// Compute total number of physical boundary faces (across all side sets)
//! \return Total number of physical boundary faces (across all side sets)
// *****************************************************************************
{
  return std::accumulate( m_bface.cbegin(), m_bface.cend(), std::size_t(0),
           []( std::size_t acc, const decltype(m_bface)::value_type& b )
              { return acc + b.second.size(); } );
}
