//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Wed 13 May 2015 09:59:41 PM MDT
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

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Performer_SDAG_CODE

  private:
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor >;

  public:
    //! Constructor
    explicit Performer( CProxy_Conductor& hostproxy,
                        LinSysMergerProxy& lsmproxy );

    //! Migrate constructor
    Performer( CkMigrateMessage* ) {}

    //! Receive matrix row contribution from fellow Performer chares
    void add(
      int id,
      const std::map< std::size_t, std::map< std::size_t, tk::real > >& rows );

  private:
    std::size_t m_id;                   //!< Charm++ array id (Base::thisIndex)
    CProxy_Conductor m_hostproxy;       //!< Host proxy
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger proxy
    std::vector< std::size_t > m_point; //!< Global ids of nodes owned
    std::map< std::size_t, std::vector< std::size_t > > m_import;
    std::map< std::size_t, std::vector< std::size_t > > m_toimport;
    std::map< std::size_t, std::map< std::size_t, tk::real > > m_lhs;

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
    initIds( const std::vector< std::size_t >& gelem ) const;

    //! Initialize points surrounding points
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
    psup() const;

    //! Assign local ids to global ids
    std::map< std::size_t, std::size_t >
    assignLid( const std::vector< std::size_t >& gid ) const;

    //! Find local for global node id
    std::size_t
    lid( const std::map< std::size_t, std::size_t >& lnode, std::size_t gid )
    const;

    //! Read coordinates of mesh nodes given
    std::array< std::vector< tk::real >, 3 >
    initCoords( const std::vector< std::size_t >& gnode ) const;

    //! Output chare mesh chare id field to file
    void
    writeChareId( const std::vector< std::size_t >& inpoel,
                  const std::array< std::vector< tk::real >, 3 >& coord ) const;

    //! Compute consistent mass matrix
    void
    consistentMass( const std::vector< std::size_t >& gnode,
                    const std::vector< std::size_t >& inpoel,
                    const std::array< std::vector< tk::real >, 3 >& coord );

    //! \brief Perform the necessary communication among fellow Performers to
    //!   update the chare-boundaries for matrix
    void
    update( const std::map< std::size_t, std::vector< std::size_t > >& exp );

    //! Contribute our portion of the left hand side matrix
    void contributeLhs();
};

} // inciter::

#endif // Performer_h
