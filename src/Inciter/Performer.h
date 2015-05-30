//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Thu 21 May 2015 10:42:05 AM MDT
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
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
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

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include <performer.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

//! Performer Charm++ chare used to advance the Euler equations in time
class Performer : public CBase_Performer {

  #if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #endif

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Performer_SDAG_CODE

  #if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic pop
  #endif

  private:
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor >;

  public:
    //! Constructor
    explicit Performer( CProxy_Conductor& hostproxy,
                        LinSysMergerProxy& lsmproxy );

    //! Migrate constructor
    explicit Performer( CkMigrateMessage* ) {}

    //! Receive matrix row contribution from fellow Performer chares
    void add(
      int id,
      const std::map< std::size_t, std::map< std::size_t, tk::real > >& rows );

  private:
    std::size_t m_id;                   //!< Charm++ array id (Base::thisIndex)
    CProxy_Conductor m_hostproxy;       //!< Host proxy
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger proxy
    std::vector< std::size_t > m_point; //!< Global ids of nodes owned
    //! Import map associating global mesh point ids to chares during import
    std::map< std::size_t, std::vector< std::size_t > > m_import;
    //! Import map associating global mesh point ids to chares before import
    std::map< std::size_t, std::vector< std::size_t > > m_toimport;
    //! Sparse matrix: global mesh point row and column ids, and nonzero value
    std::map< std::size_t, std::map< std::size_t, tk::real > > m_lhs;
    //! Right-hand side vector: global mesh point row ids and values
    std::map< std::size_t, tk::real > m_rhs;
    //! Time stamps
    std::vector< std::pair< std::string, tk::real > > m_timestamp;
    std::vector< tk::Timer > m_timer;   //!< Timers

    //! Find out if a point is owned
    bool own( std::size_t gid ) const {
      for (auto p : m_point) if (p == gid) return true;
      return false;
    }

    //! Find out if all chares have contributed we need to import from
    bool lhscomplete() const { return m_toimport == m_import; }

    //! Initialize import map
    void initImports();

    //! Initialize local->global, global->local node ids, element connectivity
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
    initIds( const std::vector< std::size_t >& gelem );

    //! Assign local ids to global ids
    std::map< std::size_t, std::size_t >
    assignLid( const std::vector< std::size_t >& gid ) const;

    //! Find local for global node id
    std::size_t
    lid( const std::map< std::size_t, std::size_t >& lnode, std::size_t gid )
    const;

    //! Read coordinates of mesh nodes given
    std::array< std::vector< tk::real >, 3 >
    initCoords( const std::vector< std::size_t >& gnode );

    //! Output chare mesh chare id field to file
    void
    writeChareId( const std::vector< std::size_t >& inpoel,
                  const std::array< std::vector< tk::real >, 3 >& coord );

    //! Compute left-hand side matrix of PDE
    void
    lhs( const std::vector< std::size_t >& gnode,
         const std::vector< std::size_t >& inpoel,
         const std::array< std::vector< tk::real >, 3 >& coord );

    //! \brief Perform the necessary communication among fellow Performers to
    //!   update the chare-boundaries for left-hand side matrix of PDE
    void
    commLhs( const std::map< std::size_t, std::vector< std::size_t > >& exp );

    //! Contribute our portion of the left-hand side matrix
    void contributeLhs();

    //! Compute righ-hand side vector of PDE
    void
    rhs( const std::vector< std::size_t >& gnode,
         const std::vector< std::size_t >& inpoel,
         const std::array< std::vector< tk::real >, 3 >& coord );
};

} // inciter::

#endif // Performer_h
