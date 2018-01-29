// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
*/
// *****************************************************************************
#ifndef FaceData_h
#define FaceData_h

#include "Types.h"
#include "Keywords.h"
#include "Fields.h"
#include "PUPUtil.h"
#include "UnsMesh.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "facedata.decl.h"

namespace tk {
  class ExodusIIMeshWriter;
  class RootMeshWriter;
}

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief FaceData class holding face-connectivity data useful
//!   to DG discretization
class FaceData
{
  public:
    //! Constructor
    explicit
      FaceData(
        const std::vector< std::size_t >& conn,
        std::size_t nbfac_complete,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::vector< std::size_t >& triinpoel_complete );

    /** @name Accessors
      * */
    ///@{
    const std::vector< std::size_t >& Gid() const { return m_gid; }
    std::vector< std::size_t >& Gid() { return m_gid; }

    const std::unordered_map< std::size_t, std::size_t >& Lid() const
    { return m_lid; }
    std::unordered_map< std::size_t, std::size_t >& Lid() { return m_lid; }

    const std::vector< std::size_t >& Inpoel() const { return m_inpoel; }
    std::vector< std::size_t >& Inpoel() { return m_inpoel; }

    const std::map< int, std::vector< std::size_t > >& Bface() const
    { return m_bface; }
    std::size_t Nbfac() const { return m_nbfac; }
    const std::vector< int >& Esuel() const { return m_esuel; }
    std::size_t Ntfac() const { return m_ntfac; }
    const std::vector< std::size_t >& Inpofa() const { return m_inpofa; }
    const std::vector< std::size_t >& Belem() const { return m_belem; }
    const std::vector< int >& Esuf() const { return m_esuf; }
    //@}

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      FaceData::pup(p);
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = std::get< 0 >( m_el );
        m_gid = std::get< 1 >( m_el );
        m_lid = std::get< 2 >( m_el );
      }
      p | m_bface;
      p | m_triinpoel;
      p | m_nbfac;
      p | m_esuel;
      p | m_ntfac;
      p | m_inpofa;
      p | m_belem;
      p | m_esuf;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i FaceData object reference
    friend void operator|( PUP::er& p, FaceData& i ) { i.pup(p); }
    //@}

  private:
    //! \brief Elements of the mesh chunk we operate on
    //! \details Initialized by the constructor. The first vector is the element
    //!   connectivity (local IDs), the second vector is the global node IDs of
    //!   owned elements, while the third one is a map of global->local node
    //!   IDs.
    std::tuple< std::vector< std::size_t >,
                std::vector< std::size_t >,
                std::unordered_map< std::size_t, std::size_t > > m_el;
    //! Alias to element connectivity in m_el
    std::vector< std::size_t > m_inpoel = std::get< 0 >( m_el );
    //! Alias to global node IDs of owned elements in m_el
    std::vector< std::size_t > m_gid = std::get< 1 >( m_el );
    //! \brief Alias to local node ids associated to the global ones of owned
    //!    elements in m_el
    std::unordered_map< std::size_t, std::size_t > m_lid = std::get< 2 >( m_el );
    //! Boundary faces side-set information
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary face-node connectivity
    std::vector< std::size_t > m_triinpoel;
    //! Number of boundary faces
    std::size_t m_nbfac;
    //! Elements surrounding elements
    std::vector< int > m_esuel;
    //! Rotal number of faces
    std::size_t m_ntfac;
    //! Face-node connectivity
    std::vector< std::size_t > m_inpofa;
    //! Boundary element vector
    std::vector< std::size_t > m_belem;
    //! Element surrounding faces
    std::vector< int > m_esuf;
};

} // inciter::

#endif // Discretization_h
