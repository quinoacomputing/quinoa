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
  std::size_t nbfac,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel ) :
  m_nbfac( nbfac ),
  m_bface( bface ),
  m_triinpoel( triinpoel )
// *****************************************************************************
//  Constructor
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//!   mesh chunk we operate on
//! \param[in] nbfac Total number of boundary-faces (triangles) in this mesh
//!   chunk
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

    for (std::size_t f=0; f<m_nbfac; ++f)
    {
      std::cout << "triinp  (" << f << "/" << m_nbfac << ") : ";
      for (std::size_t i=0; i<3; ++i)
      {
        std::cout << m_triinpoel[3*f+i] << ", ";
      }
      std::cout << "\n";
    }

    // Mapping m_triinpoel from global renumbered ids to local renumbered ids
    for (auto& i : m_triinpoel) i = tk::cref_find(lid,i);

    m_ntfac = tk::genNtfac( 4, m_nbfac, m_esuel );
    m_inpofa = tk::genInpofaTet( m_ntfac, m_nbfac, inpoel, m_triinpoel,
                                 m_esuel );
    m_belem =  tk::genBelemTet( m_nbfac, m_inpofa, esup );
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
