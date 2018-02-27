// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
*/
// *****************************************************************************
#ifndef FaceData_h
#define FaceData_h

#include "PUPUtil.h"

namespace inciter {

//! \brief FaceData class holding face-connectivity data useful
//!   to DG discretization
class FaceData
{
  public:
    //! Empty constructor for Charm++
    explicit FaceData() {}

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
