// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Transporter drives the time integration of transport equations
  \details   Transporter drives the time integration of transport equations.
*/
// *****************************************************************************
#ifndef Transporter_h
#define Transporter_h

#include <map>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "Timer.h"
#include "Types.h"
#include "InciterPrint.h"
#include "Partitioner.h"
#include "Progress.h"
#include "Scheme.h"
#include "ContainerUtil.h"

namespace inciter {

//! Transporter drives the time integration of transport equations
class Transporter : public CBase_Transporter {

  private:
    //! Indices for progress report on mesh preparation
    enum ProgMesh{ PART=0, DIST, REFINE, BND, COMM, MASK, REORD };
    //! Indices for progress report on workers preparation
    enum ProgWork{ CREATE=0, BNDFACE, COMFAC, GHOST, ADJ };

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
      #pragma clang diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
      #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__INTEL_COMPILER)
      #pragma warning( push )
      #pragma warning( disable: 1478 )
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    Transporter_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit Transporter();

    //! Reduction target: the mesh has been read from file on all PEs
    void load( std::size_t nelem, std::size_t nnode );

    //! \brief Reduction target: all Solver (PEs) have computed the number of
    //!   chares they will recieve contributions from during linear solution
    void partition();

    //! Reduction target: all PEs have distrbuted their mesh after partitioning
    void distributed();

    //! Reduction target: all PEs have created the mesh refiners
    void refinserted( int error );

    //! Reduction target: all Discretization chares have been inserted
    void discinserted();

    //! Reduction target: all Discretization constructors have been called
    void disccreated();

    //! \brief Reduction target: all worker (derived discretization) chares have
    //!   been inserted
    void workinserted();

    //! \brief Reduction target: all mesh refiner chares have setup their
    //!   boundary edges
    void edges();

    //! \brief Reduction target: all mesh refiner chares have received a round
    //!   of edges, and ran their compatibility algorithm
    void compatibility( int modified );

    //! \brief Reduction target: all mesh refiner chares have performed a step
    //!   of matching chare-boundary edges
    void matched( std::size_t nextra, std::size_t nedge, std::size_t initial );

    //! Compute surface integral across the whole problem and perform leak-test
    void bndint( tk::real sx, tk::real sy, tk::real sz );

    //! Reduction target: all PEs have optionally refined their mesh
    void refined( std::size_t nelem, std::size_t npoin );

    //! Reduction target: all Sorter chares have queried their boundary nodes
    void queried();
    //! \brief Reduction target: all Sorter chares have responded with their
    //!   boundary nodes
    void responded();

    //! \brief Reduction target: all Discretization chares have resized their
    //!    own data after mesh refinement
    void discresized();

    //! \brief Reduction target: all worker chares have resized their own data
    //!   after mesh refinement
    void workresized();

    //! Non-reduction target for receiving progress report on partitioning mesh
    void pepartitioned() { m_progMesh.inc< PART >(); }
    //! Non-reduction target for receiving progress report on distributing mesh
    void pedistributed() { m_progMesh.inc< DIST >(); }
    //! Non-reduction target for receiving progress report on finding bnd nodes
    void chbnd() { m_progMesh.inc< BND >(); }
    //! Non-reduction target for receiving progress report on node ID comm map
    void chcomm() { m_progMesh.inc< COMM >(); }
    //! Non-reduction target for receiving progress report on node ID mask
    void chmask() { m_progMesh.inc< MASK >(); }
    //! Non-reduction target for receiving progress report on reordering mesh
    void chreordered() { m_progMesh.inc< REORD >(); }

    //! Non-reduction target for receiving progress report on creating workers
    void chcreated() { m_progWork.inc< CREATE >(); }
    //! Non-reduction target for receiving progress report on finding bnd faces
    void chbndface() { m_progWork.inc< BNDFACE >(); }
    //! Non-reduction target for receiving progress report on face communication
    void chcomfac() { m_progWork.inc< COMFAC >(); }
    //! Non-reduction target for receiving progress report on sending ghost data
    void chghost() { m_progWork.inc< GHOST >(); }
    //! Non-reduction target for receiving progress report on face adjacency
    void chadj() { m_progWork.inc< ADJ >(); }

    //! Reduction target indicating that the communication maps have been setup
    void comfinal( int initial );

    //! Reduction target summing total mesh volume
    void totalvol( tk::real v, tk::real initial );

    //! \brief Reduction target indicating that all workers have finished
    //!   computing/receiving their part of the nodal volumes
    void vol();

    //! \brief Reduction target yielding the minimum mesh statistics across
    //!   all workers
    void minstat( tk::real d0, tk::real d1, tk::real d2 );

    //! \brief Reduction target yielding the maximum mesh statistics across
    //!   all workers
    void maxstat( tk::real d0, tk::real d1, tk::real d2 );

    //! \brief Reduction target yielding the sum of mesh statistics across
    //!   all workers
    void sumstat( tk::real d0, tk::real d1,
                  tk::real d2, tk::real d3,
                  tk::real d4, tk::real d5 );

    //! \brief Reduction target yielding PDF of mesh statistics across all
    //!    workers
    void pdfstat( CkReductionMsg* msg );

    //! \brief Reduction target optionally collecting diagnostics, e.g.,
    //!   residuals, from all  worker chares
    void diagnostics( CkReductionMsg* msg );

    //! Reduction target computing minimum of dt
    void advance( tk::real dt );

    //! Normal finish of time stepping
    void finish();

  private:
    InciterPrint m_print;                //!< Pretty printer
    int m_nchare;                        //!< Number of worker chares
    std::size_t m_ncit;                  //!< Number of mesh ref corr iter
    std::size_t m_nt0refit;              //!< Number of (t<0) mesh ref iters
    std::size_t m_ndtrefit;              //!< Number of (t>0) mesh ref iters
    Scheme m_scheme;                     //!< Discretization scheme
    CProxy_Partitioner m_partitioner;    //!< Partitioner nodegroup proxy
    CProxy_Refiner m_refiner;            //!< Mesh refiner array proxy
    tk::CProxy_MeshWriter m_meshwriter;  //!< Mesh writer nodegroup proxy
    CProxy_Sorter m_sorter;              //!< Mesh sorter array proxy
    std::size_t m_nelem;                 //!< Number of mesh elements
    std::size_t m_npoin_larger;          //!< Total number mesh points
     //! Total mesh volume
    tk::real m_V;
    //! Minimum mesh statistics
    std::array< tk::real, 3 > m_minstat;
    //! Maximum mesh statistics
    std::array< tk::real, 3 > m_maxstat;
    //! Average mesh statistics
    std::array< tk::real, 3 > m_avgstat;
    //! Timer tags
    enum class TimerTag { MESH_READ=0 };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;
    //! Progress object for preparing mesh
    tk::Progress< 7 > m_progMesh;
    //! Progress object for preparing workers
    tk::Progress< 5 > m_progWork;

    //! Create mesh partitioner and boundary condition object group
    void createPartitioner();

    //! Configure and write diagnostics file header
    void diagHeader();

    //! Echo diagnostics on mesh statistics
    void stat();

    //! \brief All Discretization and worker chares have resized their own data
    //!   after mesh refinement
    void resized();

    //! Query variable names for all equation systems to be integrated
    //! \param[in] eq Equation system whose variable names to query
    //! \param[in,out] var Vector of strings to which we append the variable
    //!   names for this equation. We append as many strings as many scalar
    //!   variables are in the equation system given by eq.
    template< class Eq >
    void varnames( const Eq& eq, std::vector< std::string >& var ) {
      auto o = eq.names();
      var.insert( end(var), begin(o), end(o) );
    }

    //! Verify boundary condition (BC) side sets used exist in mesh file
    //! \details This function does two things: (1) it verifies that the side
    //!   sets to which boundary conditions (BC) are assigned by the user in the
    //!   input file all exist among the side sets read from the input mesh
    //!   file and errors out if at least one does not, and (2) it matches the
    //!   side set ids at which the user has configured BCs to side set ids read
    //!   from the mesh file and removes those face and node lists associated
    //!   to side sets that the user did not set BCs on (as they will not need
    //!   processing further since they will not be used).
    //! \tparam Eq Equation type, e.g., CG, DG, used to solve PDEs
    //! \param[in] pde List of PDE system solved
    //! \param[in,out] bnd Node or face lists mapped to side set ids
    template< class Eq >
    void matchBCs( const std::vector< Eq >& pde,
                   std::map< int, std::vector< std::size_t > >& bnd )
    {
      // Query side set ids at which BCs assigned for all PDEs
      std::unordered_set< int > userbc;
      for (const auto& eq : pde) eq.side( userbc );
      // Find user-configured side set ids among side sets read from mesh file
      std::unordered_set< int > sidesets_as_bc;
      for (auto i : userbc) {   // for all side sets at which BCs are assigned
        if (bnd.find(i) != end(bnd))  // user BC found among side sets in file
          sidesets_as_bc.insert( i );  // store side set id configured as BC
        else {
          Throw( "Boundary conditions specified on side set " +
            std::to_string(i) + " which does not exist in mesh file" );
        }
      }
      // Remove sidesets not configured as BCs (will not process those further)
      using Pair = std::pair< const int, std::vector< std::size_t > >;
      tk::erase_if( bnd, [&]( Pair& item ) {
        return sidesets_as_bc.find( item.first ) == end(sidesets_as_bc);
      });
      // Warn on no BCs
      if (bnd.empty())
        m_print << "\n>>> WARNING: No boundary conditions set\n\n";
    }
};

} // inciter::

#endif // Transporter_h
