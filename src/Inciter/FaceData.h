// *****************************************************************************
/*!
  \file      src/Inciter/FaceData.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \details   Face-data used only in discontinuous Galerkin discretization scheme
*/
// *****************************************************************************
#ifndef FaceData_h
#define FaceData_h

#include <vector>
#include <tuple>
#include <map>
#include <unordered_map>

#include "Types.h"
#include "PUPUtil.h"
#include "ContainerUtil.h"

namespace inciter {

//! Data associated to a tetrahedron cell id used to comunicate across faces
using GhostData =
  std::unordered_map< std::size_t, // tet id
                      std::tuple<
                        // 3 node ids for potentially multiple faces
                        std::vector< std::size_t >,
                        // elem geometry, see tk::genGeoElemTet()
                        std::vector< tk::real >,
                        // coordinates of vertex of tet that is not on face
                        std::array< tk::real, 3 >,
                        // relative position of above vertex in inpoel
                        std::size_t,
                        // inpoel of said tet
                        std::array< std::size_t, 4 > > >;


//! FaceData class holding face-connectivity data useful for DG discretization
class FaceData {

  public:
    //! Empty constructor for Charm++
    explicit FaceData() {}

    //! Constructor
    explicit
    FaceData( const std::vector< std::size_t >& ginpoel,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::map< int, std::vector< std::size_t > >& bnode,
              std::vector< std::size_t >& triinpoel );

    /** @name Accessors
      * */
    ///@{
    const std::map< int, std::vector< std::size_t > >& Bface() const
    { return m_bface; }
    const std::map< int, std::vector< std::size_t > >& Bnode() const
    { return m_bnode; }
    std::map< int, std::vector< std::size_t > >& Bnode() { return m_bnode; }
    const std::vector< std::size_t >& Triinpoel() const { return m_triinpoel; }
    std::size_t Nbfac() const { return tk::sumvalsize( m_bface ); }
    const std::vector< int >& Esuel() const { return m_esuel; }
    std::size_t Ntfac() const { return m_ntfac; }
    const std::vector< std::size_t >& Inpofa() const { return m_inpofa; }
    std::vector< std::size_t >& Inpofa() { return m_inpofa; }
    const std::vector< std::size_t >& Belem() const { return m_belem; }
    const std::vector< int >& Esuf() const { return m_esuf; }
    std::vector< int >& Esuf() { return m_esuf; }
    //@}

    /** @name Charm++ pack/unpack (serialization) routines
      * */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | m_bface;
      p | m_bnode;
      p | m_triinpoel;
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
    //! Boundary nodes side-set information
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Triangle face connecitivity
    std::vector< std::size_t > m_triinpoel;
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
