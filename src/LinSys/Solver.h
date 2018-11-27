// *****************************************************************************
/*!
  \file      src/LinSys/Solver.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ linear system merger nodegroup to solve a linear system
  \details   Charm++ linear system merger nodegroup used to collect and
    assemble the left hand side matrix (lhs), the right hand side (rhs) vector,
    and the solution (unknown) vector from individual worker
    chares. Beside collection and assembly, the system is also solved. The
    solution is outsourced to hypre, an MPI-only library. Once the solution is
    available, the individual worker chares are updated with the new solution.

    This class assembles and solves two linear systems, whose rhs vectors may
    change during time stepping. One of the two linear systems is called a
    high-order and the other one is the low-order linear system.

    Characteristics of the low-order linear system: (1) the left hand side
    matrix is diagonal, (2) the right hand side vector is a combination of the
    high-order right hand side vector and another vector, assembled separately
    (and overlapped with that of the high-order system). This dual system
    solution is done here as required by the flux-corrected transport algorithm
    used for transport equations in inciter, and the low order system is
    derived from the high order system using continuous Galerkin node-centered
    scheme.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef Solver_h
#define Solver_h

#include <vector>
#include <map>
#include <unordered_map>

#include "Types.h"
#include "Tags.h"
#include "Fields.h"
#include "TaggedTuple.h"
#include "HypreMatrix.h"
#include "HypreVector.h"
#include "HypreSolver.h"
#include "Callback.h"

#include "NoWarning/solver.decl.h"

namespace tk {

//! Linear system merger and solver Charm++ chare nodegroup class
//! \details Instantiations of Solver comprise a processor aware Charm++
//!   chare nodegroup. When instantiated, a new object is created on each
//!   compute node and not more (as opposed to individual chares or chare
//!   array object elements). The nodegroup's elements are used to collect
//!   information from all chare objects that happen to be on a given compute
//!   node. See also the Charm++ interface file solver.ci.
class Solver : public CBase_Solver {

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #pragma GCC diagnostic ignored "-Wunused-parameter"
  #elif defined(__INTEL_COMPILER)
    #pragma warning( push )
    #pragma warning( disable: 1478 )
  #endif
  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Solver_SDAG_CODE
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  public:
    //! Constructor
    Solver( const SolverCallback& cb, std::size_t n );

    //! Set number of worker chares expected to contribute on this compute node
    void nchare( int n );

    //!  Receive lower and upper global node IDs from chares
    void chbounds( std::size_t lower, std::size_t upper );

    //!  Communicate lower and upper bounds across all compute nodes
    void nodebounds( int n, std::size_t lower, std::size_t upper );

    //! Chares contribute their ids and callbacks
    void charecom( int fromch, const MatCGCallback& cb );

    //! Chares contribute their global row ids for establishing communications
    void charerow( int fromch, const std::vector< std::size_t >& row );

    //! Prepare for next step
    void next();

    //! Receive global row ids from fellow nodegroup branches
    void addrow( int fromch, const std::set< std::size_t >& row );

    //! Chares contribute their solution nonzero values
    void charesol( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& solution );

    //! Receive solution vector nonzeros from fellow nodegroup branches
    void addsol( int fromch,
                 const std::map< std::size_t,
                                 std::vector< tk::real > >& solution );

    //! Chares contribute their matrix nonzero values
    void charelhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const std::pair< std::vector< std::size_t >,
                                    std::vector< std::size_t > >& psup,
                   const tk::Fields& lhsd,
                   const tk::Fields& lhso );

    //! Receive matrix nonzeros from fellow nodegroup branches
    void addlhs( int fromch,
                 const std::map< std::size_t,
                                 std::map< std::size_t,
                                           std::vector< tk::real > > >& l );

    //! Chares contribute their rhs nonzero values
    void charerhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& r );

    //! Receive right-hand side vector nonzeros from fellow nodegroup branches
    void addrhs( int fromch,
                 const std::map< std::size_t, std::vector< tk::real > >& r );

    //! Chares contribute to the rhs of the low-order linear system
    void charelowrhs( int fromch,
                      const std::vector< std::size_t >& gid,
                      const Fields& lowrhs );

    //! Receive low-order rhs vector nonzeros from fellow node group branches
    void addlowrhs( int fromch,
                    const std::map< std::size_t,
                                    std::vector< tk::real > >& lowrhs );

    //! Chares contribute to the lhs of the low-order linear system
    void charelowlhs( int fromch,
                      const std::vector< std::size_t >& gid,
                      const Fields& lowlhs );

    //! \brief Receive lhs vector to the low-order system from fellow
    //!   nodegroup branches
    void addlowlhs( int fromch,
                    const std::map< std::size_t,
                                    std::vector< tk::real > >& lowlhs );

    //! All communications have been establised among compute nodes
    void comfinal();

    //! Chares query Dirichlet boundary conditions
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same compute node.
    const std::unordered_map< std::size_t,
            std::vector< std::pair< bool, tk::real > > >&
      dirbc() { return m_bc; }

    //! \brief Chares contribute their global row ids and associated Dirichlet
    //!   boundary condition values at which they set BCs
    void charebc( const std::unordered_map< std::size_t,
                          std::vector< std::pair< bool, tk::real > > >& bc );

    //! Reduction target collecting the final aggregated BC node list map
    void addbc( CkReductionMsg* msg );

    //! \brief Chares contribute their numerical and analytical solutions
    //!   nonzero values for computing diagnostics
    void charediag( int fromch,
                    uint64_t it,
                    tk::real t,
                    tk::real dt,
                    const std::vector< std::size_t >& gid,
                    const Fields& u,
                    const Fields& a,
                    const std::vector< tk::real >& v );

    //! \brief Receive numerical and analytical solution vector nonzeros from
    //!   fellow nodegroup branches for computing diagnostics
    void adddiag( int fromch,
                  std::map< std::size_t,
                    std::vector< std::vector< tk::real > > >& solution );

  private:
    SolverCallback m_cb;    //!< Charm++ associated to compile-time tags
    std::size_t m_ncomp;    //!< Number of scalar components per unknown
    std::size_t m_nchare;   //!< Total number of worker chares
    std::size_t m_mynchare; //!< Number of chares contributing to my node
    std::size_t m_nbounds;  //!< Number of chares contributed bounds to my node
    std::size_t m_ncomm;    //!< Number of chares finished commaps on my node
    std::size_t m_nchbc;    //!< Number of chares we received bcs from
    std::size_t m_lower;    //!< Lower index of the global rows on my node
    std::size_t m_upper;    //!< Upper index of the global rows on my node
    uint64_t m_it;          //!< Iteration count (original in Discretization)
    tk::real m_t;           //!< Physical time (original in Discretization)
    tk::real m_dt;          //!< Time step size (original in Discretization)
    bool m_initial;         //!< True in the first time step, false after
    //! Chare id and callbacks to entry methods of all worker chares
    std::map< int, tk::MatCGCallback > m_worker;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the global row ids
    std::vector< std::vector< std::size_t > > m_rowimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the solution/unknown vector
    std::vector< std::vector< std::size_t > > m_solimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the left-hand side matrix
    std::vector< std::vector< std::size_t > > m_lhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the righ-hand side vector
    std::vector< std::vector< std::size_t > > m_rhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the low-order rhs vector
    std::vector< std::vector< std::size_t > > m_lowrhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the low-order lhs vector
    std::vector< std::vector< std::size_t > > m_lowlhsimport;
    //! Part of global row indices owned by my node for each contributing chare
    std::vector< std::set< std::size_t > > m_row;
    //! \brief Part of unknown/solution vector owned by my node
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row IDs
    std::vector< std::vector< tk::real > > m_sol;
    //! \brief Part of left-hand side matrix owned by my node
    //! \details Nonzero values (for each scalar equation solved) associated to
    //!   global mesh point row and column IDs.
    std::map< std::size_t,
              std::map< std::size_t, std::vector< tk::real > > > m_lhs;
    //! \brief Part of right-hand side vector owned by my node
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids
    std::map< std::size_t, std::vector< tk::real > > m_rhs;
    //! \brief Part of low-order right-hand side vector owned by my node
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the rhs
    //!   terms to be combined with the rhs to produce the low-order rhs.
    std::map< std::size_t, std::vector< tk::real > > m_lowrhs;
    //! Part of the low-order system left-hand side vector owned by my node
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the nonzero values
    //!   of the low-order system lhs "matrix" solution.
    std::map< std::size_t, std::vector< tk::real > > m_lowlhs;
    tk::hypre::HypreVector m_x; //!< Hypre vector to store the solution/unknowns
    tk::hypre::HypreMatrix m_A; //!< Hypre matrix to store the left-hand side
    tk::hypre::HypreVector m_b; //!< Hypre vector to store the right-hand side
    //! Hypre solver
    tk::hypre::HypreSolver m_solver;
    //! Row indices for my node
    std::vector< int > m_hypreRows;
    //! Number of matrix columns/rows on my node
    std::vector< int > m_hypreNcols;
    //! Matrix column indices for rows on my node
    std::vector< int > m_hypreCols;
    //! Matrix nonzero values for my node
    std::vector< tk::real > m_hypreMat;
    //! RHS vector nonzero values for my node
    std::vector< tk::real > m_hypreRhs;
    //! Solution vector nonzero values for my node
    std::vector< tk::real > m_hypreSol;
    //! Compute nodes associated to lower and upper global row indices
    //! \details These are the divisions at which the linear system is divided
    //!   along compute node boundaries.
    std::map< std::pair< std::size_t, std::size_t >, int > m_div;
    //! Compute nodes associated to global mesh point indices
    //! \details This is used to cache the compute node associated to mesh nodes
    //!   communicated, so that a quicker-than-linear-cost search can be used to
    //!   find the compute node for a communicated mesh node after the node is
    //!   in the cache.
    std::map< std::size_t, int > m_node;
    //! \brief Values (for each scalar equation solved) of Dirichlet boundary
    //!   conditions assigned to global node IDs we set
    //! \details The map key is the global mesh node/row ID, the value is a
    //!   vector of pairs in which the bool (.first) indicates whether the
    //!   boundary condition value (.second) is set at the given node. The size
    //!   of the vectors is the number of PDEs integrated times the number of
    //!   scalar components in all PDEs. This BC map stores all row IDs at which
    //!   Dirichlet boundary conditions are prescribed, i.e., including boundary
    //!   conditions set across all compute nodes, not just the ones need to be
    //!   set on this compute node.
    std::unordered_map< std::size_t,
                        std::vector< std::pair< bool, tk::real > > > m_bc;
    //! Matrix non-zero coefficients for Dirichlet boundary conditions
    std::unordered_map< std::size_t,
                        std::map< std::size_t, std::vector<tk::real> > > m_bca;

    //! Return compute node id for global mesh row id
    int node( std::size_t gid );

    //! Serialize solution vector into a Charm++ message, ready for a CkCallback
    std::pair< int, std::unique_ptr<char[]> >
    serializeSol( const std::vector< std::size_t >& gid,
                  const std::vector< tk::real >& u ) const;

    //! Check if our portion of the solution vector values is complete
    bool solcomplete() const { return m_solimport == m_rowimport; }

    //! Check if our portion of the matrix values is complete
    bool lhscomplete() const { return m_lhsimport == m_rowimport; }

    //! Check if our portion of the right-hand side vector values is complete
    bool rhscomplete() const { return m_rhsimport == m_rowimport; }

    //! Check if our portion of the low-order rhs vector values is complete
    bool lowrhscomplete() const { return m_lowrhsimport == m_rowimport; }

    //! Check if our portion of the low-order lhs vector values is complete
    bool lowlhscomplete() const { return m_lowlhsimport == m_rowimport; }

    //! Build Hypre data for our portion of the global row ids
    void hyprerow();

    //! Set boundary conditions on the left-hand side matrix
    void lhsbc();

    //! Set Dirichlet boundary conditions on the right-hand side vector
    void rhsbc();

    //! Build Hypre data for our portion of the solution vector
    void hypresol();

    //! Build Hypre data for our portion of the matrix
    void hyprelhs();

    //! Build Hypre data for our portion of the right-hand side vector
    void hyprerhs();

    //! Set our portion of values of the distributed solution vector
    void sol();

    //! Set our portion of values of the distributed matrix
    void lhs();

    //! Set our portion of values of the distributed right-hand side vector
    void rhs();

    //! Assemble distributed solution vector
    void assemblesol();

    //! Assemble distributed matrix
    void assemblelhs();

    //! Assemble distributed right-hand side vector
    void assemblerhs();

    //! Update solution vector in workers contributing on this compute node
    void updateSol();

    //! Solve hyigh-order linear system
    void solve();

    //! \brief Update low-order solution vector in workers contributing on this
    //!   compute node
    void updateLowSol();

    //! Solve low-order linear system
    void lowsolve();

    //! Update diagnostics vector
    void updateDiag( std::size_t row,
                     std::vector< tk::real >&& u,
                     std::vector< tk::real >&& a,
                     tk::real v );

    //! Compute diagnostics (residuals) and contribute them back to host
    void diagnostics();
};

} // tk::

#endif // Solver_h
