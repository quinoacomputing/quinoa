// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
  \see       FaceData.h for more info.
*/
// *****************************************************************************

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
  std::size_t nbfac_complete,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel_complete ) : m_bface( bface )
// *****************************************************************************
//  Constructor
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//!   mesh chunk we operate on
//! \param[in] nbfac_complete Total number of boundary-faces (triangles) in 
//!   the entire mesh
//! \param[in] bface Map of boundary-face lists mapped to corresponding 
//!   side set ids
//! \param[in] triinpoel_complete Interconnectivity of points and 
//!   boundary-face in the entire mesh
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
    m_nbfac = tk::genNbfacTet( nbfac_complete, inpoel, triinpoel_complete,
                               m_triinpoel );
    m_esuel = tk::genEsuelTet( inpoel,tk::genEsup(inpoel,4) );
    m_ntfac = tk::genNtfac( 4, m_nbfac, m_esuel );
    m_inpofa = tk::genInpofaTet( m_ntfac, m_nbfac, inpoel, m_triinpoel,
                                 m_esuel );
    m_belem =  tk::genBelemTet( m_nbfac, m_inpofa, tk::genEsup(inpoel,4) );
    m_esuf = tk::genEsuf( 4, m_ntfac, m_nbfac, m_belem, m_esuel );

    Assert( m_belem.size() == m_nbfac,
            "Number of boundary-elements and number of boundary-faces unequal" );
  }
}

// void
// FaceData::test()
// // *****************************************************************************
// // test member function
// // *****************************************************************************
// {
// }
