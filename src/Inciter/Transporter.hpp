// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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

#include "Timer.hpp"
#include "Types.hpp"
#include "InciterPrint.hpp"
#include "Partitioner.hpp"
#include "Progress.hpp"
#include "Scheme.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

//! Indices for progress report on mesh preparation
enum ProgMesh{ PART=0, DIST, REFINE, BND, COMM, MASK, REORD };
//! Prefixes for progress report on mesh preparation
static const std::array< std::string, 7 >
  ProgMeshPrefix = {{ "p", "d", "r", "b", "c", "m", "r" }},
  ProgMeshLegend = {{ "partition", "distribute", "refine", "bnd", "comm",
                      "mask", "reorder" }};

//! Indices for progress report on workers preparation
enum ProgWork{ CREATE=0, BNDFACE, COMFAC, GHOST, ADJ };
//! Prefixes for progress report on workers preparation
static const std::array< std::string, 5 >
  ProgWorkPrefix = {{ "c", "b", "f", "g", "a" }},
  ProgWorkLegend = {{ "create", "bndface", "comfac", "ghost", "adj" }};

//! Transporter drives the time integration of transport equations
class Transporter : public CBase_Transporter {

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

    //! Migrate constructor: returning from a checkpoint
    explicit Transporter( CkMigrateMessage* m );

    //! Reduction target: the mesh has been read from file on all PEs
    void load( std::size_t meshid, std::size_t nelem );

    //! Reduction target: a mesh has been partitioned
    void partitioned( std::size_t meshid );

    //! \brief Reduction target: all Solver (PEs) have computed the number of
    //!   chares they will recieve contributions from during linear solution
    void partition();

    //! \brief Reduction target: all compute nodes have distrbuted their mesh
    //!   after partitioning
    void distributed( std::size_t meshid );

    //! Reduction target: all Refiner chares have queried their boundary edges
    void queriedRef( std::size_t meshid );
    //! Reduction target: all Refiner chares have setup their boundary edges
    void respondedRef( std::size_t meshid );

    //! Reduction target: all compute nodes have created the mesh refiners
    void refinserted( std::size_t meshid, std::size_t error );

    //! Reduction target: all Discretization chares have been inserted
    void discinserted( std::size_t meshid );

    //! Reduction target: all Discretization constructors have been called
    void disccreated( std::size_t summeshid, std::size_t npoin );

    //! \brief Reduction target: all worker (derived discretization) chares have
    //!   been inserted
    void workinserted( std::size_t meshid );

    //! \brief Reduction target: all Refiner chares have received a round
    //!   of edges, and have run their compatibility algorithm
    void compatibility( std::size_t meshid );

    //! \brief Reduction target: all Refiner chares have matched/corrected
    //!   the tagging of chare-boundary edges, all chares are ready to perform
    //!   refinement.
    void matched( std::size_t summeshid, std::size_t nextra, std::size_t nref,
                  std::size_t nderef, std::size_t sumrefmode );

    //! Compute surface integral across the whole problem and perform leak-test
    void bndint( tk::real sx, tk::real sy, tk::real sz, tk::real cb,
                 tk::real summeshid );

    //! Reduction target: all chares have refined their mesh
    void refined( std::size_t summeshid, std::size_t nelem, std::size_t npoin );

    //! \brief Reduction target: all worker chares have resized their own data
    //!   after AMR or ALE
    void resized( std::size_t meshid );

    //! Reduction target: all worker chares have generated their own esup
    void startEsup( std::size_t meshid );

    //! Reduction target: all Sorter chares have queried their boundary edges
    void queried( std::size_t meshid );
    //! Reduction target: all Sorter chares have setup their boundary edges
    void responded( std::size_t meshid );

    //! Non-reduction target for receiving progress report on partitioning mesh
    void pepartitioned() { m_progMesh.inc< PART >( printer() ); }
    //! Non-reduction target for receiving progress report on distributing mesh
    void pedistributed() { m_progMesh.inc< DIST >( printer() ); }
    //! Non-reduction target for receiving progress report on finding bnd nodes
    void chbnd() { m_progMesh.inc< BND >( printer() ); }
    //! Non-reduction target for receiving progress report on node ID comm map
    void chcomm() { m_progMesh.inc< COMM >( printer() ); }
    //! Non-reduction target for receiving progress report on node ID mask
    void chmask() { m_progMesh.inc< MASK >( printer() ); }
    //! Non-reduction target for receiving progress report on reordering mesh
    void chreordered() { m_progMesh.inc< REORD >( printer() ); }

    //! Non-reduction target for receiving progress report on creating workers
    void chcreated() { m_progWork.inc< CREATE >( printer() ); }
    //! Non-reduction target for receiving progress report on finding bnd faces
    void chbndface() { m_progWork.inc< BNDFACE >( printer() ); }
    //! Non-reduction target for receiving progress report on face communication
    void chcomfac() { m_progWork.inc< COMFAC >( printer() ); }
    //! Non-reduction target for receiving progress report on sending ghost data
    void chghost() { m_progWork.inc< GHOST >( printer() ); }
    //! Non-reduction target for receiving progress report on face adjacency
    void chadj() { m_progWork.inc< ADJ >( printer() ); }

    // Reduction target indicating all "ghosts" insertions are done
    void doneInsertingGhosts(std::size_t meshid);

    //! Reduction target indicating that the communication maps have been setup
    void comfinal( std::size_t initial, std::size_t summeshid );

    //! Reduction target summing total mesh volume
    void totalvol( tk::real v, tk::real initial, tk::real summeshid );

    //! \brief Reduction target yielding the minimum mesh statistics across
    //!   all workers
    void minstat( tk::real d0, tk::real d1, tk::real d2, tk::real rmeshid );

    //! \brief Reduction target yielding the maximum mesh statistics across
    //!   all workers
    void maxstat( tk::real d0, tk::real d1, tk::real d2, tk::real rmeshid );

    //! \brief Reduction target yielding the sum of mesh statistics across
    //!   all workers
    void sumstat( tk::real d0, tk::real d1,
                  tk::real d2, tk::real d3,
                  tk::real d4, tk::real d5,
                  tk::real summeshid );

    //! \brief Reduction target yielding PDF of mesh statistics across all
    //!    workers
    void pdfstat( CkReductionMsg* msg );

    //! Reduction target computing total volume of IC box
    void boxvol( tk::real v, tk::real summeshid );

    //! \brief Reduction target optionally collecting diagnostics, e.g.,
    //!   residuals, from all  worker chares
    void diagnostics( CkReductionMsg* msg );

    //! Resume execution from checkpoint/restart files
    void resume();

    //! Save checkpoint/restart files
    void checkpoint( std::size_t finished, std::size_t meshid );

    //! Normal finish of time stepping
    void finish( std::size_t meshid = 0 );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ mainchare, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_input;
      p | m_nchare;
      p | m_meshid;
      p | m_nload;
      p | m_npart;
      p | m_nstat;
      p | m_ndisc;
      p | m_nchk;
      p | m_ncom;
      p | m_ncit;
      p | m_nt0refit;
      p | m_ndtrefit;
      p | m_noutrefit;
      p | m_noutderefit;
      p | m_scheme;
      p | m_partitioner;
      p | m_refiner;
      p | m_meshwriter;
      p | m_sorter;
      p | m_nelem;
      p | m_npoin;
      // returning from checkpoint
      if (p.isUnpacking()) m_finished.resize( m_nchare.size(), 0 );
      p | m_meshvol;
      p | m_minstat;
      p | m_maxstat;
      p | m_avgstat;
      p | m_timer;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] t Transporter object reference
    friend void operator|( PUP::er& p, Transporter& t ) { t.pup(p); }
    //@}

  private:
    //! List of mesh files to be used for potentially multiple solvers
    std::vector< std::string > m_input;
    //! Number of worker chares (one per mesh)
    std::vector< int > m_nchare;
    //! Sum of mesh ids (across all chares, key) for each meshid (value)
    std::unordered_map< std::size_t, std::size_t > m_meshid;
    //! Number of mesh ref corr iter (one per mesh)
    std::vector< std::size_t > m_ncit;
    //! Number of meshes loaded
    std::size_t m_nload;
    //! Number of meshes partitioned
    std::size_t m_npart;
    //! Number of mesh statistics computed
    std::size_t m_nstat;
    //! Number of Discretization arrays created
    std::size_t m_ndisc;
    //! Number of worker arrays checkpointed
    std::size_t m_nchk;
    //! Number of worker arrays have finished setting up their comm maps
    std::size_t m_ncom;
    //! Number of t0ref mesh ref iters (one per mesh)
    std::vector< std::size_t > m_nt0refit;
    //! Number of dtref mesh ref iters (one per mesh)
    std::vector< std::size_t > m_ndtrefit;
    //! Number of outref mesh ref iters (one per mesh)
    std::vector< std::size_t > m_noutrefit;
    //! Number of outderef mesh ref iters (one per mesh)
    std::vector< std::size_t > m_noutderefit;
    //! Discretization scheme (one per mesh)
    std::vector< Scheme > m_scheme;
    //! Partitioner nodegroup proxies (one per mesh)
    std::vector< CProxy_Partitioner > m_partitioner;
    //! Mesh refiner array proxies (one per mesh)
    std::vector< CProxy_Refiner > m_refiner;
    //! Mesh writer nodegroup proxies (one per mesh)
    std::vector< tk::CProxy_MeshWriter > m_meshwriter;
    //! Mesh sorter array proxy (one per mesh)
    std::vector< CProxy_Sorter > m_sorter;
    //!< Number of mesh elements (per mesh)
    std::vector< std::size_t > m_nelem;
    //!< Number of mesh points (per mesh)
    std::vector< std::size_t > m_npoin;
    //!< Nonzero if finished with timestepping (one per mesh)
    std::vector< std::size_t > m_finished;
    //! Total mesh volume (one per mesh)
    std::vector< tk::real > m_meshvol;
    //! Minimum mesh statistics (one per mesh)
    std::vector< std::array< tk::real, 3 > > m_minstat;
    //! Maximum mesh statistics (one per mesh)
    std::vector< std::array< tk::real, 3 > > m_maxstat;
    //! Average mesh statistics (one per mesh)
    std::vector< std::array< tk::real, 3 > > m_avgstat;
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

    //! Echo configuration to screen
    void info( const InciterPrint& print );

    //! Print out time integration header to screen
    void inthead( const InciterPrint& print );

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

    //! Create pretty printer specialized to Inciter
    //! \return Pretty printer
    InciterPrint printer() const {
      const auto& def =
        g_inputdeck_defaults.get< tag::cmd, tag::io, tag::screen >();
      auto nrestart = g_inputdeck.get< tag::cmd, tag::io, tag::nrestart >();
      return InciterPrint(
        g_inputdeck.get< tag::cmd >().logname( def, nrestart ),
        g_inputdeck.get< tag::cmd, tag::verbose >() ? std::cout : std::clog,
        std::ios_base::app );
    }

    //! Function object for querying the side set ids the user configured
    //! \details Used to query and collect the side set ids the user has
    //!   configured for all PDE types querying all BC types. Used on a
    //!   Carteisan product of 2 type lists: PDE types and BC types.
    struct UserBC {
      const ctr::InputDeck& inputdeck;
      std::unordered_set< int >& userbc;
      explicit UserBC( const ctr::InputDeck& i, std::unordered_set< int >& u )
        : inputdeck(i), userbc(u) {}
      template< typename U > void operator()( brigand::type_<U> ) {
        using tag::param;
        using eq = typename brigand::front< U >;
        using bc = typename brigand::back< U >;
        for (const auto& s : inputdeck.get< param, eq, tag::bc, bc >())
          for (const auto& i : s) userbc.insert( std::stoi(i) );
      }
    };

    //! Function object for querying the side set ids for time dependent BCs
    //! \details Used to query and collect the side set ids the user has
    //!   configured for all PDE types querying time dependent BCs. Used on
    //!   PDE type list.
    struct UserTimedepBC {
      const ctr::InputDeck& inputdeck;
      std::unordered_set< int >& userbc;
      explicit UserTimedepBC( const ctr::InputDeck& i,
        std::unordered_set< int >& u )
        : inputdeck(i), userbc(u) {}
      template< typename eq > void operator()( brigand::type_<eq> ) {
        using tag::param;
        for (const auto& sys : inputdeck.get< param, eq, tag::bctimedep >()) {
          for (const auto& b : sys) {
            for (auto i : b.template get< tag::sideset >())
              userbc.insert( std::stoi(i) );
          }
        }
      }
    };

    //! Verify boundary condition (BC) side sets used exist in mesh file
    bool matchBCs( std::map< int, std::vector< std::size_t > >& bnd );

    //! Print out mesh statistics
    void meshstat( const std::string& header ) const;

    //! Generate list of input mesh filenames configured by the user
    std::vector< std::string > input();

    //! Decide if ALE will need a linear solver
    bool need_linearsolver() const;
};

} // inciter::

#endif // Transporter_h
