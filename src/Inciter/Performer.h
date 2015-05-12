//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Tue 12 May 2015 07:38:59 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances the Euler equations
  \details   Performer advances the Euler equations. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances the Euler equations in time.
*/
//******************************************************************************
#ifndef Performer_h
#define Performer_h

#include <iostream>     // NOT REALLY NEEDED

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <inciter.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <DerivedData.h>
#include <Inciter/InputDeck/InputDeck.h>
#include <ExodusIIMeshReader.h>
#include <ExodusIIMeshWriter.h>
#include <LinSysMerger.h>
#include <performer.decl.h>

namespace inciter {

//! Performer Charm++ chare used to advance the Euler equations in time
class Performer : public CBase_Performer {

  private:
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor >;

  public:
    //! Constructor
    explicit Performer( CProxy_Conductor& hostproxy,
                        LinSysMergerProxy& lsmproxy );

    //! Migrate constructor
    Performer( CkMigrateMessage* ) {}

  private:
    std::size_t m_id;                   //!< Charm++ array id (Base::thisIndex)
    CProxy_Conductor m_hostproxy;       //!< Host proxy
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger proxy
    std::vector< std::size_t > m_point; //!< Global ids of nodes owned
    //! Export map
    std::map< std::size_t, std::vector< std::size_t > > m_export;
    //! Points surrounding points
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;

    //! Initialize local->global, global->local node ids, element connectivity
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
    initIds( const std::vector< std::size_t >& gelem );

    //! Initialize points surrounding points
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > psup();

    //! Assign local ids to global ids
    std::map< std::size_t, std::size_t >
    assignLid( const std::vector< std::size_t >& gid );

    //! Find local for global node id
    std::size_t
    lid( const std::map< std::size_t, std::size_t >& lnode, std::size_t gid );

    //! Read coordinates of mesh nodes given
    std::array< std::vector< tk::real >, 3 >
    initCoords( const std::vector< std::size_t >& gnode );

    //! Output chare mesh chare id field to file
    void writeChareId( const std::vector< std::size_t >& inpoel,
                       const std::array< std::vector< tk::real >, 3 >& coord );
};

} // inciter::

#endif // Performer_h
