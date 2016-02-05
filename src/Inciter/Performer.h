//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Fri 05 Feb 2016 06:25:46 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances a PDE
  \details   Performer advances a PDE. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a PDE in time.
*/
//******************************************************************************
#ifndef Performer_h
#define Performer_h

#include <array>
#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>
#include <cstring>
#include <cmath>
#include <unordered_map>

#include "Types.h"
#include "MeshNodes.h"
#include "Inciter/InputDeck/InputDeck.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "conductor.decl.h"
#include "performer.decl.h"
#include "linsysmerger.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace tk { class ExodusIIMeshWriter; }

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Performer Charm++ chare used to advance a PDE in time
class Performer : public CBase_Performer {

  private:
    using ConductorProxy = CProxy_Conductor;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor,
                                                       CProxy_Performer >;
  public:
    //! Constructor
    explicit
      Performer( const CProxy_Conductor& conductor,
                 const LinSysMergerProxy& lsm,
                 const std::vector< std::size_t >& conn,
                 const std::unordered_map< std::size_t, std::size_t >& cid );

    //! Migrate constructor
    explicit Performer( CkMigrateMessage* ) {}

    //! Initialize mesh IDs, element connectivity, coordinates
    void setup();

    //! Initialize communication and mesh data
    void init( tk::real dt );

    //! Update solution vector
    void updateSolution( const std::vector< std::size_t >& gid,
                         const std::vector< tk::real >& sol );

    //! Advance equations to next stage in multi-stage time stepping
    void advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t );

  private:
    uint64_t m_it;                      //!< Iteration count
    uint64_t m_itf;                     //!< Field output iteration count
    tk::real m_t;                       //!< Physical time
    uint8_t m_stage;                    //!< Stage in multi-stage time stepping
    std::size_t m_nsol;                 //!< Counter for solution nodes updated
    ConductorProxy m_conductor;         //!< Conductor proxy
    LinSysMergerProxy m_linsysmerger;   //!< Linear system merger proxy
    //! \brief Map associating old node IDs (as in file) to new node IDs (as in
    //!   producing contiguous-row-id linear system contributions)
    std::unordered_map< std::size_t, std::size_t > m_cid;
    std::vector< std::size_t > m_inpoel;//!< Element connectivity (local IDs)
    std::vector< std::size_t > m_gid;   //!< Global node ids of owned elements
    //!< Local node ids associated to the global ones of owned elements
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Mesh point coordinates
    std::array< std::vector< tk::real >, 3 > m_coord;
    //! Unknown/solution vector: global mesh point row ids and values
    std::vector< tk::real > m_u, m_uf, m_un;

    //! Send off global row IDs to linear system merger, setup global->local IDs
    void setupIds();

    //! Read coordinates of mesh nodes given
    void readCoords();

    //! Set initial conditions
    void ic();

    //! Compute left-hand side matrix of PDE
    void lhs();

    //! Compute righ-hand side vector of PDE
    void rhs( tk::real mult,
              tk::real dt,
              const std::vector< tk::real >& sol,
              std::vector< tk::real >& rhs );

    //! Output chare mesh to file
    void writeMesh();

    //! Output chare mesh chare id field to file
    void writeChareId( const tk::ExodusIIMeshWriter& ew, uint64_t it ) const;

    //! Output solution to file
    void writeSolution( const tk::ExodusIIMeshWriter& ew,
                        uint64_t it,
                        int varid,
                        const std::vector< tk::real >& u ) const;

    //! Output mesh-based fields metadata to file
    void writeMeta() const;

    //! Output mesh-based fields to file
    void writeFields( tk::real time );

    //! Compute initial conditions for dispersion in simple shear flow
    tk::real ansol_shear( std::size_t i, tk::real t ) const {
      const tk::real X0 = 7200.0;        // x position of source
      const tk::real t0 = g_inputdeck.get< tag::discr, tag::t0 >(); // initial t
      const tk::real U0 = 0.5;           // velocity in x direction
      const tk::real LAMBDA = 5.0e-4;    // amount of shear, (lambda = du/dy)
      const tk::real D = 10.0;           // scalar diffusivity
      const tk::real M = 4.0*M_PI*t0*std::sqrt( 1.0 + LAMBDA*LAMBDA*t0*t0/12.0 );
      const auto& x = m_coord[0];
      const auto& y = m_coord[1];
      tk::real a = x[i] - X0 - U0*t - 0.5*LAMBDA*y[i]*t;
      tk::real b = 1.0 + LAMBDA*LAMBDA*t*t/12.0;
      return M * exp( -a*a/(4.0*M_PI*D*t*b) - y[i]*y[i]/(4.0*D*t) )
               / ( 4.0*M_PI*t*std::sqrt(b) );
    };

    //! Compute initial conditions representing a joint Gaussian
    tk::real ansol_gauss( std::size_t i ) const {
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
