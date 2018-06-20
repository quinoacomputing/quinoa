// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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
#include "VectorReducer.h"
#include "Progress.h"
#include "Scheme.h"

namespace inciter {

//! Transporter drives the time integration of transport equations
class Transporter : public CBase_Transporter {

  private:
    //! Indices for progress report on mesh preparation
    enum ProgMesh{ PART=0, DIST, REFINE, BND, QUERY, MASK, REORD, BOUND };

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
    void load( uint64_t nelem );

    //! \brief Reduction target: all Solver (PEs) have computed the number of
    //!   chares they will recieve contributions from during linear solution
    void nchare();

    //! Reduction target: all Solver (PEs) have computed their row bounds
    void bounds();

    //! Reduction target: all PEs have distrbuted their mesh after partitioning
    void distributed();

    //! Reduction target: all PEs have created the mesh refiners
    void refinserted();

    //! Reduction target: all Discretization chares have been inserted
    void discinserted();

    //! Reduction target: all Discretization constructors have been called
    void disccreated();

    //! \brief Reduction target: all worker (drevied discretization) chares have
    //!   been inserted
    void workinserted();

    //! \brief Reduction target: all mesh refiner chares have distributed their
    //!   newly added node IDs that are shared among chares
    void matched( std::size_t extra );

    //! Reduction target: all PEs have optionally refined their mesh
    void refined( std::size_t nelem, std::size_t npoin );

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished flattening its global mesh node IDs and they are ready
    //!   for computing the communication maps required for node ID reordering
    void flattened();

    //! Non-reduction target for receiving progress report on partitioning mesh
    void pepartitioned() { m_progMesh.inc< PART >(); }
    //! Non-reduction target for receiving progress report on distributing mesh
    void pedistributed() { m_progMesh.inc< DIST >(); }
    //! Non-reduction target for receiving progress report on mesh refinement
    void chrefined() { m_progMesh.inc< REFINE >(); }
    //! Non-reduction target for receiving progress report on flattening mesh
    void chbnd() { m_progMesh.inc< BND >(); }
    //! Non-reduction target for receiving progress report on node ID query
    void chquery() { m_progMesh.inc< QUERY >(); }
    //! Non-reduction target for receiving progress report on node ID mask
    void chmask() { m_progMesh.inc< MASK >(); }
    //! Non-reduction target for receiving progress report on reordering mesh
    void chreordered() { m_progMesh.inc< REORD >(); }
    //! Non-reduction target for receiving progress report on computing bounds
    void chbounds() { m_progMesh.inc< BOUND >(); }

    //! \brief Reduction target indicating that the communication has been
    //!    established among PEs
    void comfinal();

    //! Reduction target summing total mesh volume
    void totalvol( tk::real v );

    //! \brief Reduction target indicating that all workers have finished
    //!   computing/receiving their part of the nodal volumes
    void vol();

    //! \brief Reduction target yielding the minimum mesh statistics across
    //!   all workers
    void minstat( tk::real d0, tk::real d1 );

    //! \brief Reduction target yielding the maximum mesh statistics across
    //!   all workers
    void maxstat( tk::real d0, tk::real d1 );

    //! \brief Reduction target yielding the sum of mesh statistics across
    //!   all workers
    void sumstat( tk::real d0, tk::real d1, tk::real d2, tk::real d3 );

    //! \brief Reduction target yielding PDF of mesh statistics across all
    //!    workers
    void pdfstat( CkReductionMsg* msg );

    //! \brief Reduction target optionally collecting diagnostics, e.g.,
    //!   residuals, from all  worker chares
    void diagnostics( CkReductionMsg* msg );

    //! Start time stepping
    void start();

    //! Reset linear solver for next time step
    void next() { m_solver.next(); }

    //! Normal finish of time stepping
    void finish();

  private:
    InciterPrint m_print;                //!< Pretty printer
    int m_nchare;                        //!< Number of worker chares
    uint64_t m_chunksize;                //!< Number of elements per PE
    uint64_t m_remainder;                //!< Number elements added to last PE
    tk::CProxy_Solver m_solver;          //!< Linear system solver group proxy
    Scheme m_scheme;                     //!< Discretization scheme
    CProxy_Partitioner m_partitioner;    //!< Partitioner group proxy
    CProxy_Refiner m_refiner;            //!< Mesh refiner array proxy
    CProxy_Sorter m_sorter;              //!< Mesh sorter array proxy
    std::size_t m_nelem;                 //!< Number mesh elements
    std::size_t m_npoin;                 //!< Number mesh points
    //! Average communication cost of merging the linear system
    tk::real m_avcost;
     //! Total mesh volume
    tk::real m_V;
    //! Minimum mesh statistics
    std::array< tk::real, 2 > m_minstat;
    //! Maximum mesh statistics
    std::array< tk::real, 2 > m_maxstat;
    //! Average mesh statistics
    std::array< tk::real, 2 > m_avgstat;
    //! Timer tags
    enum class TimerTag { MESH_READ=0 };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;
    //! \brief Aggregate 'old' (as in file) node ID list at which Solver
    //!   sets boundary conditions, see also Partitioner.h
    std::vector< std::size_t > m_linsysbc;
    //! Progress object for preparing mesh
    tk::Progress< 8 > m_progMesh;

    //! Create linear solver group
    void createSolver();

    //! Create mesh partitioner and boundary condition object group
    void createPartitioner();

    //! Start partitioning the mesh
    void partition();

    //! Configure and write diagnostics file header
    void diagHeader();

    //! Echo diagnostics on mesh statistics
    void stat();

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

    //! Verify boundarty condition (BC) side sets used exist in mesh file
    //! \details This function verifies that the side sets to which boundary
    //!   conditions B(C) are assigned by the user in the input file all exist
    //!   in the mesh file and warns if at least one does not.
    //! \tparam Eq Equation type, e.g., CG, DG, we are using to solver PDEs
    //! \param[in] pde List of PDE system solved
    //! \param[in,out] er ExodusII mesh reader object
    //! \note Failure here is not an error, bot only a wanring, i.e., the user
    //!   must check the screen output for this warning. Would this be more user
    //!   friendly to make it an error, i.e., abort with an error message?
    //! \note If the input vector is empty, no checking is done.
    template< class Eq >
    void verifyBCsExist( const std::vector< Eq >& pde,
                         tk::ExodusIIMeshReader& er )
    {
      // Query BC assigned to all side sets for all PDEs
      std::unordered_set< int > conf;
      for (const auto& eq : pde) eq.side( conf );
      // Read in side sets associated to mesh node IDs from file
      auto sidenodes = er.readSidesets();
      for (auto i : conf) {
        if (sidenodes.find(i) == end(sidenodes)) {
          Throw( "WARNING: Boundary conditions specified on side set " +
                 std::to_string(i) + " which does not exist in mesh file" );
          break;
        }
      }
    }
};

} // inciter::

#endif // Transporter_h
