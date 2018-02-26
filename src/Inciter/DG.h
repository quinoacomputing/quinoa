// *****************************************************************************
/*!
  \file      src/Inciter/DG.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     DG advances a system of PDEs with the discontinuous Galerkin scheme
  \details   DG advances a system of partial differential equations (PDEs) using
    discontinuous Galerkin (DG) finite element (FE) spatial discretization (on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping.

    There are a potentially large number of DG Charm++ chares created by
    Transporter. Each DG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication.
*/
// *****************************************************************************
#ifndef DG_h
#define DG_h

#include <vector>

#include "DerivedData.h"
#include "FaceData.h"

#include "NoWarning/dg.decl.h"

namespace inciter {

//! DG Charm++ chare array used to advance PDEs in time with DG+RK
class DG : public CBase_DG {

  public:
    //! Constructor
    explicit DG( const CProxy_Discretization& disc,
                 const tk::CProxy_Solver& solver,
                 const FaceData& fd );

    //! Migrate constructor
    explicit DG( CkMigrateMessage* ) {}

    //! Configure Charm++ reduction types for concatenating BC nodelists
    static void registerReducers();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup( tk::real v );

    //! Compute time step size
    void dt();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Evaluate whether to continue with next step
    void eval();

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_DG::pup(p);
      p | m_itf;
      p | m_disc;
      p | m_fd;
      p | m_u;
      p | m_un;
      p | m_vol;
      p | m_geoFace;
      p | m_geoElem;
      p | m_ax;
      p | m_ay;
      p | m_az;
      p | m_lhs;
      p | m_rhs;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i DG object reference
    friend void operator|( PUP::er& p, DG& i ) { i.pup(p); }
    //@}

  private:
    //! Field output iteration count
    uint64_t m_itf;
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Face data
    FaceData m_fd;
    //! Vector of unknown/solution average over each mesh element
    tk::Fields m_u;
    //! Vector of unknown at previous time-step
    tk::Fields m_un;
    //! Total mesh volume
    tk::real m_vol;
    //! Face geometry
    tk::Fields m_geoFace;
    //! Element geometry
    tk::Fields m_geoElem;
    //! advection velocity x, y and z components
    tk::real m_ax;
    tk::real m_ay;
    tk::real m_az;
    //! Left-hand side mass-matrix which is a diagonal matrix
    std::vector< tk::real > m_lhs;
    //! Vector of right-hand side
    tk::Fields m_rhs;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Compute left hand side
    void lhs();

    //! Compute right hand side
    void rhs();

    //! Upwind fluxes
    std::vector< tk::real > upwindFlux( std::vector< tk::real > ul,
                                        std::vector< tk::real > ur,
                                        std::array< tk::real, 3 > fn );

    //! Time stepping
    void solve( tk::real dt );

    //! Prepare for next step
    void next();
 
    //! Output mesh and particle fields to files
    void out();

    //! Compute diagnostics, e.g., residuals
    bool diagnostics();

    //! Output mesh-based fields to file
    void writeFields( tk::real time );
};

} // inciter::

#endif // DG_h
