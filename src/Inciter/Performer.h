//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Sat 30 May 2015 11:41:21 AM MDT
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

#include <array>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <utility>
#include <vector>
#include <cstring>

#include "Timer.h"
#include "Types.h"
#include "LinSysMerger.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "inciter.decl.h"
#include "performer.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace tk { class ExodusIIMeshWriter; }

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
    void addLhs(
      int id,
      const std::map< std::size_t, std::map< std::size_t, tk::real > >& rows );

    //! Receive right-hand side vector contribution from fellow Performer chares
    void addRhs( int id, const std::map< std::size_t, tk::real >& rows );

  private:
    std::size_t m_id;                   //!< Charm++ array id (Base::thisIndex)
    CProxy_Conductor m_hostproxy;       //!< Host proxy
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger proxy
    int m_it;                           //!< Iteration count
    tk::real m_t;                       //!< Physical time
    std::vector< std::size_t > m_point; //!< Global ids of nodes owned
    std::size_t m_nelem;                //!< Number of owned elements
    //! Import map associating global mesh point ids to chares during lhs import
    std::map< std::size_t, std::vector< std::size_t > > m_lhsimport;
    //! Import map associating global mesh point ids to chares during rhs import
    std::map< std::size_t, std::vector< std::size_t > > m_rhsimport;
    //! Import map associating global mesh point ids to chares before import
    std::map< std::size_t, std::vector< std::size_t > > m_toimport;
    //! Sparse matrix: global mesh point row and column ids, and nonzero value
    std::map< std::size_t, std::map< std::size_t, tk::real > > m_lhs;
    //! Right-hand side vector: global mesh point row ids and values
    std::map< std::size_t, tk::real > m_rhs;
    //! Unknown vector: global mesh point row ids and values
    std::map< std::size_t, tk::real > m_x;
    //! Time stamps
    std::vector< std::pair< std::string, tk::real > > m_timestamp;
    enum class TimerTag { LHS, RHS, SOL };     //!< Timer labels
    std::map< TimerTag, tk::Timer > m_timer;   //!< Timers

    //! Find out if a point is owned
    bool own( std::size_t gid ) const {
      for (auto p : m_point) if (p == gid) return true;
      return false;
    }

    //! Find out if all chares have contributed we need to import lhs parts from
    bool lhscomplete() const { return m_toimport == m_lhsimport; }

    //! Find out if all chares have contributed we need to import rhs parts from
    bool rhscomplete() const { return m_toimport == m_rhsimport; }

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

    //! Set initial conditions
    void
    ic( const std::vector< std::size_t >& gnode,
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

    //! \brief Perform the necessary communication among fellow Performers to
    //!   update the chare-boundaries for right-hand side vector of PDE
    void
    commRhs( const std::map< std::size_t, std::vector< std::size_t > >& exp );

    //! Contribute our portion of the left-hand side matrix
    void contributeLhs();

    //! Contribute our portion of the right-hand side vector
    void contributeRhs();

    //! Compute righ-hand side vector of PDE
    void
    rhs( const std::vector< std::size_t >& gnode,
         const std::vector< std::size_t >& inpoel,
         const std::array< std::vector< tk::real >, 3 >& coord );

    //! Output chare mesh to file
    void
    writeMesh( const std::vector< std::size_t >& inpoel,
               const std::array< std::vector< tk::real >, 3 >& coord );

    //! Output chare mesh chare id field to file
    void
    writeChareId( const tk::ExodusIIMeshWriter& ew ) const;

    //! Output solution to file
    void writeSolution( const tk::ExodusIIMeshWriter& ew ) const;

    //! Output mesh-based fields metadata to file
    void writeMeta() const;

    //! Output mesh-based fields to file
    void writeFields();
};

} // inciter::

#endif // Performer_h
