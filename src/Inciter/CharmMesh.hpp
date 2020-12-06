// *****************************************************************************
/*!
  \file      src/Inciter/CharmMesh.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Chare class declaration for holding part of a mesh
  \details   Chare class declaration for holding part of a mesh, for
             interoperating with Charm++ libraries.
*/
// *****************************************************************************
#ifndef CharmMesh_h
#define CharmMesh_h

#include "UnsMesh.hpp"

#include "NoWarning/charm++.hpp"
#include "NoWarning/charmmesh.decl.h"

namespace inciter {

//! CharmMesh chare array holding part of a mesh
class CharmMesh : public CBase_CharmMesh {

  public:
    //! Constructor
    explicit CharmMesh( const std::vector< std::size_t >& inpoel,
                        const tk::UnsMesh::Coords& coord,
                        const tk::Fields& u,
                        int nchare );

    //! Pass source mesh to transfer library
    void transferSource();

    //! Pass destination mesh to transfer library
    void transferDest();

    //! Mesh transfer complete
    void solutionFound();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_inpoel;
      p | m_coord;
      p | m_u;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] m CharmMesh object reference
    friend void operator|( PUP::er& p, CharmMesh& m ) { m.pup(p); }
    //@}

  private:
    //! Mesh connectivity graph
    std::vector< std::size_t > m_inpoel;
    //! Mesh node coordinates
    tk::UnsMesh::Coords m_coord;
    //! Numerical solution to transfer
    tk::Fields m_u;
};

} // inciter::

#endif // CharmMesh_h
