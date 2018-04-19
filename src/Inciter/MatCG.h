// *****************************************************************************
/*!
  \file      src/Inciter/MatCG.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     MatCG for a PDE system with continuous Galerkin with a matrix
  \details   MatCG advances a system of partial differential equations (PDEs)
    using continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a time
    stepping scheme that is equivalent to the Lax-Wendroff (LW) scheme within
    the unstructured-mesh FE context and treats discontinuities with
    flux-corrected transport (FCT). The left-hand side (consistent-mass) matrix
    is stored in a compressed sparse row (CSR) storage and thus this scheme uses
    a matrix-based linear solver.

    There are a potentially large number of MatCG Charm++ chares created by
    Transporter. Each MatCG gets a chunk of the full load (part of the mesh)
    and does the same: initializes and advances a number of PDE systems in time.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/matcg.ci.
*/
// *****************************************************************************
#ifndef MatCG_h
#define MatCG_h

#include <vector>
#include <map>
#include <unordered_set>

#include "QuinoaConfig.h"
#include "Types.h"
#include "Fields.h"
#include "DerivedData.h"
#include "VectorReducer.h"
#include "FluxCorrector.h"
#include "NodeDiagnostics.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "FaceData.h"

#include "NoWarning/matcg.decl.h"

namespace tk {
  class ExodusIIMeshWriter;
  class RootMeshWriter;
}

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! MatCG Charm++ chare array used to advance PDEs in time with MatCG+LW+FCT
class MatCG : public CBase_MatCG {

  public:
    //! Constructor
    explicit MatCG( const CProxy_Discretization& disc,
                    const tk::CProxy_Solver& solver,
                    const FaceData& fd );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit MatCG( CkMigrateMessage* ) : m_diag( *Disc() ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ custom reduction types initiated from this chare array
    static void registerReducers();

    //! Setup: query boundary conditions, output mesh, etc.
    void setup( tk::real v );

    //! Compute time step size
    void dt();

    //! Advance equations to next time step
    void advance( tk::real newdt );

    //! Update high order solution vector
    void updateSol( //solMsg* m );
                    const std::vector< std::size_t >& gid,
                    const std::vector< tk::real >& sol );

    //! Update low order solution vector
    void updateLowSol( //solMsg* m );
                       const std::vector< std::size_t >& gid,
                       const std::vector< tk::real >& sol );

    //! Prepare for next step
    void next( const tk::Fields& a );

    //! Evaluate whether to continue with next step
    void eval();

    //! Verify that solution does not change at Dirichlet boundary conditions
    bool correctBC( const tk::Fields& a );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_MatCG::pup(p);
      p | m_itf;
      p | m_nhsol;
      p | m_nlsol;
      p | m_disc;
      p | m_solver;
      p | m_u;
      p | m_ul;
      p | m_du;
      p | m_dul;
      p | m_ue;
      p | m_lhsd;
      p | m_lhso;
      p | m_vol;
      p | m_diag;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i MatCG object reference
    friend void operator|( PUP::er& p, MatCG& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Field output iteration count
    uint64_t m_itf;
    //! Counter for high order solution nodes updated
    std::size_t m_nhsol;
    //! Counter for low order solution nodes updated
    std::size_t m_nlsol;
    //! Discretization proxy
    CProxy_Discretization m_disc;
    //! Linear system merger and solver proxy
    tk::CProxy_Solver m_solver;
    //! \brief Map associating local node IDs to side set IDs for all side sets
    //!   read from mesh file (not only those the user sets BCs on)
    std::map< int, std::vector< std::size_t > > m_side;
    //! Unknown/solution vector at mesh nodes
    tk::Fields m_u;
    //! Unknown/solution vector at mesh nodes (low orderd)
    tk::Fields m_ul;
    //! Unknown/solution vector increment (high order)
    tk::Fields m_du;
    //! Unknown/solution vector increment (low order)
    tk::Fields m_dul;
    //! Unknown/solution vector at mesh cells
    tk::Fields m_ue;
    //! Sparse matrix sotring the diagonals and off-diagonals of nonzeros
    tk::Fields m_lhsd, m_lhso;
    //! Total mesh volume
    tk::real m_vol;
    //! Diagnostics object
    NodeDiagnostics m_diag;

    //! Access bound Discretization class pointer
    Discretization* Disc() const {
      Assert( m_disc[ thisIndex ].ckLocal() != nullptr, "ckLocal() null" );
      return m_disc[ thisIndex ].ckLocal();
    }

    //! Output mesh and particle fields to files
    void out();

    //! Output mesh-based fields to file
    void writeFields( tk::real time );

    //! \brief Extract node IDs from side set node lists and match to
    //    user-specified boundary conditions
    void bc();

    //! Compute left-hand side of transport equations
    void lhs();

    //! Compute righ-hand side vector of transport equations
    void rhs();
};

} // inciter::

#endif // MatCG_h
