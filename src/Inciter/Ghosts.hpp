// *****************************************************************************
/*!
  \file      src/Inciter/Ghosts.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Declarations file for generating ghost data structures
  \details   Declarations file for asynchronous distributed
             ghost data structures using Charm++.

    There are a potentially large number of Ghosts Charm++ chares.
    Each Ghosts chare gets a chunk of the full load, due to partiting the mesh.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef Ghosts_h
#define Ghosts_h

#include "ConjugateGradients.hpp"
#include "Inciter/Options/MeshVelocitySmoother.hpp"
#include "Options/UserTable.hpp"
#include "Fields.hpp"
#include "Table.hpp"
#include "FaceData.hpp"

#include "NoWarning/ghosts.decl.h"

namespace inciter {

//! Ghosts Charm++ chare array used to determine ghost data structures
class Ghosts : public CBase_Ghosts {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
      #pragma clang diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
      #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__INTEL_COMPILER)
      #pragma warning( push )
      #pragma warning( disable: 1478 )
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    Ghosts_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit
    Ghosts( const std::map< int, std::vector< std::size_t > >& bface,
      const std::vector< std::size_t >& triinpoel );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Ghosts( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    /** @name Pack/unpack (Charm++ serialization) routines */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_disc;
      p | m_fd;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] a Ghosts object reference
    friend void operator|( PUP::er& p, Ghosts& a ) { a.pup(p); }
    ///@}

  private:
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Face data
    FaceData m_fd;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }
};

} // inciter::

#endif // Ghosts_h
