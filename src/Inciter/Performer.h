//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Sun 14 Jun 2015 09:35:02 PM MDT
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

#include "Types.h"
#include "Timer.h"
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

//   #if defined(__clang__) || defined(__GNUC__)
//     #pragma GCC diagnostic push
//     #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//   #endif
// 
//   // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
//   // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
//   Performer_SDAG_CODE
// 
//   #if defined(__clang__) || defined(__GNUC__)
//     #pragma GCC diagnostic pop
//   #endif

  private:
    using LinSysMergerProxy =
      tk::CProxy_LinSysMerger< CProxy_Conductor, CProxy_Performer >;

  public:
    //! Constructor
    explicit Performer( CProxy_Conductor& hostproxy,
                        LinSysMergerProxy& lsmproxy );

    //! Migrate constructor
    explicit Performer( CkMigrateMessage* ) {}

    //! Build linear system by computing matrix, unknown and rhs vectors
    void buildSystem();

    //! Update solution vector
    void updateSolution( const std::map< std::size_t, tk::real >& sol );

//     //! Receive matrix row contribution from fellow Performer chares
//     void addLhs(
//       int id,
//       const std::map< std::size_t, std::map< std::size_t, tk::real > >& rows );

//     //! Receive right-hand side vector contribution from fellow Performer chares
//     void addRhs( int id, const std::map< std::size_t, tk::real >& rows );

  private:
    std::size_t m_id;                   //!< Charm++ array id (Base::thisIndex)
    CProxy_Conductor m_hostproxy;       //!< Host proxy
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger proxy
    int m_it;                           //!< Iteration count
    tk::real m_t;                       //!< Physical time
    std::vector< std::size_t > m_point; //!< Global ids of nodes owned
    std::vector< std::size_t > m_gid;   //!< Global node ids of owned elements
    std::vector< std::size_t > m_inpoel;//!< Owned element connectivity
    //! Mesh point coordinates
    std::array< std::vector< tk::real >, 3 > m_coord;
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
    //! Unknown/solution vector: global mesh point row ids and values
    std::map< std::size_t, tk::real > m_sol;
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
    void initIds( const std::vector< std::size_t >& gelem );

    //! Assign local ids to global ids
    std::map< std::size_t, std::size_t >
    assignLid( const std::vector< std::size_t >& gid ) const;

    //! Find local for global node id
    std::size_t
    lid( const std::map< std::size_t, std::size_t >& lnode, std::size_t gid )
    const;

    //! Read coordinates of mesh nodes given
    void initCoords();

    //! Set initial conditions
    void ic();

    //! Compute left-hand side matrix of PDE
    void lhs();

//     //! \brief Perform the necessary communication among fellow Performers to
//     //!   update the chare-boundaries for left-hand side matrix of PDE
//     void
//     commLhs( const std::map< std::size_t, std::vector< std::size_t > >& exp );

//     //! \brief Perform the necessary communication among fellow Performers to
//     //!   update the chare-boundaries for right-hand side vector of PDE
//     void
//     commRhs( const std::map< std::size_t, std::vector< std::size_t > >& exp );

//     //! Contribute our portion of the left-hand side matrix
//     void contributeLhs();

//     //! Contribute our portion of the right-hand side vector
//     void contributeRhs();

    //! Compute righ-hand side vector of PDE
    void rhs();

    //! Output chare mesh to file
    void writeMesh();

    //! Output chare mesh chare id field to file
    void
    writeChareId( const tk::ExodusIIMeshWriter& ew ) const;

    //! Output solution to file
    void writeSolution( const tk::ExodusIIMeshWriter& ew ) const;

    //! Output mesh-based fields metadata to file
    void writeMeta() const;

    //! Output mesh-based fields to file
    void writeFields();

    //! Compute initial conditions for dispersion in simple shear flow
    tk::real
    ansol_shear( std::size_t i, tk::real t ) const {
      const tk::real X0 = 0.5;           // x position of source
      const tk::real Y0 = 0.5;           // y position of source
      const tk::real U0 = 0.5;           // velocity in x direction
      const tk::real t0 = 2.4;           // initial time
      const tk::real LAMBDA = 5.0e-4;    // amount of shear, (lambda = du/dy)
      const tk::real D = 10.0;           // scalar diffusivity
      const tk::real MASS = 4 * M_PI * t0 *
                            std::sqrt( 1.0 + LAMBDA * LAMBDA * t0 * t0 / 12.0 );
      const auto& x = m_coord[0];
      const auto& y = m_coord[1];
      tk::real a = x[i] - X0 - U0*t - 0.5*LAMBDA*y[i]*t;
      tk::real b = 1.0 + LAMBDA*LAMBDA*t*t/12.0;
      return MASS * exp( -a*a/(4*D*t*b) - y[i]*y[i]/(4*D*t) )
                  / ( 4.0 * M_PI * t * std::sqrt(b) );
    };

    //! Compute initial conditions representing a joint Gaussian
    tk::real
    ansol_gauss( std::size_t i ) const {
      const tk::real X0 = 0.5, Y0 = 0.5, Z0 = 0.5;
      const tk::real v1 = 0.1, v2 = 0.05, v3 = 0.05;
      const auto& x = m_coord[0];
      const auto& y = m_coord[1];
      const auto& z = m_coord[2];
      const tk::real rx = x[i]-X0;
      const tk::real ry = y[i]-Y0;
      const tk::real rz = z[i]-Z0;
      return std::exp( -0.5 * (rx*rx/v1 + ry*ry/v2 + rz*rz/v3) ) /
             2.0 / M_PI / v1 / v2 / v3;
    }

};

} // inciter::

#endif // Performer_h
