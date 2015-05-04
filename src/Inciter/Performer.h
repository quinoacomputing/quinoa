//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Mon 04 May 2015 09:14:51 AM MDT
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

    //! Merge performer linear system contributions to PEs
    void initLinearSystem();

    //! Return our processing element number to caller's group branch
    void pe( int branch ) { m_lsmproxy[ branch ].workPe( m_id, CkMyPe() ); }

  private:
    std::size_t m_id;                   //!< Charm++ array id (Base::thisIndex)
    CProxy_Conductor m_hostproxy;       //!< Host proxy
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger proxy

    //! Global ids of nodes owned
    std::vector< std::size_t > m_point;

    //! Number of owned mesh points
    std::size_t m_nown;

    //! Communication maps
    std::map< std::size_t, std::vector< std::size_t > > m_export;
    //std::map< std::size_t, std::vector< std::size_t > > m_import;

    //! Points surrounding points
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;

    //! Initialize import map
    void initImports();

    //! Initialize local->global, global->local node ids, element connectivity
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
    initIds( const std::vector< std::size_t >& gelem );

    //! Initialize data structures derived from mesh connectivity
    void initDerivedData();

    //! Assign local ids to global ids
    std::map< std::size_t, std::size_t >
    assignLid( const std::vector< std::size_t >& gid );

    //! Find local for global node id
    std::size_t
    lid( const std::map< std::size_t, std::size_t >& lnode, std::size_t gid );

    //! Read coordinates of owned and received mesh nodes
    std::array< std::vector< tk::real >, 3 >
    initCoords( const std::vector< std::size_t >& gnode );

    //! Register ourselves with our PE's linear system merger
    void registerWithLinSysMerger();

    //! Output chare mesh and nodal chare id field to file
    void writeChareId( const std::vector< std::size_t >& inpoel,
                       const std::array< std::vector< tk::real >, 3 >& coord );
};

} // inciter::

#endif // Performer_h
