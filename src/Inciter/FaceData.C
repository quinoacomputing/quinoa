// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
  \see       FaceData.h for more info.
*/
// *****************************************************************************

#include "Tags.h"
#include "Reorder.h"
#include "Vector.h"
#include "DerivedData.h"
#include "Discretization.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Scheme.h"
#include "PDE.h"
#include "Print.h"

namespace inciter {

static CkReduction::reducerType PDFMerger;
extern std::vector< PDE > g_pdes;
extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::FaceData;

FaceData::FaceData(
  const std::vector< std::size_t >& conn,
  std::size_t nbfac_complete,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel_complete ) :
  (g_inputdeck.get< tag::selected, tag::scheme >() == ctr::SchemeType::DG) ? (
    m_el( tk::global2local( conn ) ),     // fills m_inpoel, m_gid, m_lid
    m_bface( bface ),
    m_nbfac( tk::genNbfacTet( nbfac_complete, m_inpoel, triinpoel_complete, m_triinpoel ) ),
    m_esuel( tk::genEsuelTet( m_inpoel,tk::genEsup(m_inpoel,4) ) ),
    m_ntfac( tk::genNtfac( 4, m_nbfac, m_esuel ) ),
    m_inpofa( tk::genInpofaTet( m_ntfac, m_nbfac, m_inpoel, m_triinpoel, m_esuel ) ),
    m_belem( tk::genBelemTet( m_nbfac, m_inpofa, tk::genEsup(m_inpoel,4) ) ),
    m_esuf( tk::genEsuf( 4, m_ntfac, m_nbfac, m_belem, m_esuel ) )
  ) : ( {} )
// *****************************************************************************
//  Constructor
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//!   mesh chunk we operate on
//! \details "Contiguous-row-id" here means that the numbering of the mesh nodes
//!   (which corresponds to rows in the linear system) are (approximately)
//!   contiguous (as much as this can be done with an unstructured mesh) as the
//!   problem is distirbuted across PEs, held by Solver objects. This ordering
//!   is in start contrast with "as-in-file" ordering, which is the ordering of
//!   the mesh nodes as it is stored in the file from which the mesh is read in.
//!   The as-in-file ordering is highly non-contiguous across the distributed
//!   problem.
// *****************************************************************************
{
  Assert( m_belem.size() == m_nbfac,
          "Number of boundary-elements and number of boundary-faces unequal" );
}

void
FaceData::test()
// *****************************************************************************
// test member function
// *****************************************************************************
{
}

#include "facedata.def.h"
