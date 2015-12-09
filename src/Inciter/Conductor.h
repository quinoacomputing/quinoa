//******************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Wed 09 Dec 2015 08:22:56 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Conductor drives the time integration of a PDE
  \details   Conductor drives the time integration of a PDE
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/conductor.ci.
*/
//******************************************************************************
#ifndef Conductor_h
#define Conductor_h

#include <map>
#include <vector>
#include <iosfwd>
#include <utility>

#include "Timer.h"
#include "Types.h"
#include "InciterPrint.h"
#include "Partitioner.h"
#include "Performer.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "conductor.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

//! Conductor drives the time integration of a PDE
class Conductor : public CBase_Conductor {

  public:
    //! Constructor
    explicit Conductor();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished reading their part of the computational mesh graph and
    //!   we are ready to compute the computational load
    void load( uint64_t nelem );

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished setting up the necessary data structures for
    //!   partitioning the computational mesh and we are ready for partitioning
    void partition();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished distributing the mesh element IDs after partitioning and
    //!   we ready to start reordering mesh node IDs whose first step is reading
    //!   the mesh graph connectivity (i.e., the global node IDs) owned on each
    //!   PE
    void readOwnedGraph();

    //! Reorder mesh node IDs owned on each PE
    void reorder( int pe );

    //! \brief Charm++ entry method inidicating that all Partitioner chare
    //!  groups have finished preparing the computational mesh
    void prepared(
      int pe,
      const std::unordered_map< int, std::vector< std::size_t > >& conn,
      const std::unordered_map< int,
              std::unordered_map< std::size_t, std::size_t > >& chcid );

    //! Start a round of estimation by querying the communication costs from PEs
    void query(
      const std::map< int, std::pair< std::size_t, std::size_t > >& div );

    //! \brief Receive and estimate overall communication cost of merging the
    //!   linear system
    void costed( int pe, tk::real c );

    //! \brief Reduction target indicating that all Spawner chare groups have
    //!   finished creating their Charm++ Performer worker chare array elements
    //!   (initializing their mesh element ids they will work on)
    void created();

    //! \brief Reduction target indicating that all linear system merger
    //!   branches have done their part of storing and exporting global row ids
    void rowcomplete();

    //! \brief Reduction target indicating that all Performer chares have
    //!   finished their initialization step
    void initcomplete();

    //! \brief Reduction target indicating that all Performer chares have
    //!   finished a time step and it is time to decide whether to continue
    void evaluateTime();

    //! \brief Reduction target indicating that all ...
    void advance();

    //! Normal finish of time stepping
    void finish();

    //! Collect vector of time stamps from (Performer) chare array
    //! \param[in] stamp Vector of time stamps contributed    
    void arrTimestamp(
      const std::vector< std::pair< std::string, tk::real > >& stamp )
    { timestamp( stamp, m_arrTimestamp, m_arrTimestampCnt, m_nchare ); }

    //! Collect vector of time stamps from (LinSysMerger) chare group branches
    //! \param[in] stamp Vector of time stamps contributed
    void grpTimestamp(
      const std::vector< std::pair< std::string, tk::real > >& stamp )
    { timestamp( stamp, m_grpTimestamp, m_grpTimestampCnt, CkNumPes() ); }

    //! Collect performance statistics from (Performer) chare array elements
    //! \param[in] p Vector of performance statistics contributed    
    void arrPerfstat(
      const std::vector< std::pair< std::string, tk::real > >& p )
    { perfstat( p, m_arrPerfstat, m_arrPerfstatCnt, m_nchare ); }

    //! Collect performance statistics from (LinSysMerger) chare group branches
    //! \param[in] p Vector of performance statistics contributed
    void grpPerfstat(
      const std::vector< std::pair< std::string, tk::real > >& p )
    { perfstat( p, m_grpPerfstat, m_grpPerfstatCnt, CkNumPes() ); }

  private:
    using PartitionerProxy = CProxy_Partitioner< CProxy_Conductor >;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor,
                                                       CProxy_Performer >;
    using SpawnerProxy = CProxy_Spawner< CProxy_Conductor,
                                         CProxy_Performer,
                                         LinSysMergerProxy >;

    InciterPrint m_print;               //!< Pretty printer
    int m_nchare;                       //!< Number of performer chares
    int m_prepared;                     //!< perpared() chare group counter
    std::size_t m_costed;               //!< costed() chare group counter
    int m_eval;                         //!< EvaluateTime() charge group counter
    int m_init;                         //!< initcomplete() charge group counter
    uint64_t m_it;                      //!< Iteration count
    tk::real m_t;                       //!< Physical time
    tk::real m_dt;                      //!< Physical time step size
    uint8_t m_stage;                    //!< Stage in multi-stage time stepping
    int m_arrTimestampCnt;              //!< Time stamp chare array counter
    int m_grpTimestampCnt;              //!< Time stamp chare group counter
    int m_arrPerfstatCnt;               //!< Perfstat chare array counter
    int m_grpPerfstatCnt;               //!< Perfstat chare group counter
    PartitionerProxy m_partitioner;     //!< Partitioner group
    SpawnerProxy m_spawner;             //!< Spawner group
    LinSysMergerProxy m_linsysmerger;   //!< Linear system merger chare group
    //! \brief  Reordered global mesh connectivity, i.e., node IDs, all PEs
    //!   (outer key) with their chares (inner key) contribute to in a linear
    //!   system
    std::unordered_map< int,
      std::unordered_map< int, std::vector< std::size_t > > > m_conn;
    //! \brief Map associating old node IDs (as in file) to new node IDs (as
    //!   in producing contiguous-row-id linear system contributions)
    //!   categorized by chare IDs (middle key) associated to PEs (outer key)
    std::unordered_map< int,
      std::unordered_map< int,
        std::unordered_map< std::size_t, std::size_t > > > m_chcid;
    //! \brief Lower and upper global row indices associated to a PE
    //! \details These are the divisions at which the linear system is divided
    //!   at along PE boundaries.
    std::map< int, std::pair< std::size_t, std::size_t > > m_div;
    //! Communication cost for merging the linear system associated to PE
    std::map< int, tk::real > m_cost;
    //! Time stamps merged from chare array elements
    std::map< std::string, std::vector< tk::real > > m_arrTimestamp;
    //! Time stamps merged from chare group elements
    std::map< std::string, std::vector< tk::real > > m_grpTimestamp;
    //! Performance statistics merged from chare array elements
    std::map< std::string, std::vector< tk::real > > m_arrPerfstat;
    //! Performance statistics merged from chare group elements
    std::map< std::string, std::vector< tk::real > > m_grpPerfstat;
    //! Timer labels
    enum class TimerTag { MESH,
                          CREATE,
                          SETUP,
                          INITIALIZE,
                          TIMESTEP };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;

    //! Query communication cost after mesh reordering
    void computeCost();

    //! Estimate communication cost across all PEs
    std::pair< tk::real, tk::real >
    estimate( const std::map< int, tk::real >& cost );

    //! Create linear system merger group and worker chares
    void createWorkers();

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header();

    //! Print out one-liner report on time step
    void report();

    //! Send collected timer and performance data to host
    void finalReport();

    //! Collect and compute averages of time stamps contributed by chares
    void timestamp(
      const std::vector< std::pair< std::string, tk::real > >& stamp,
      std::map< std::string, std::vector< tk::real > >& map,
      int& counter,
      int max );

    //! Collect and compute performance statistics contributed by chares
    void perfstat( const std::vector< std::pair< std::string, tk::real > >& p,
                   std::map< std::string, std::vector< tk::real > >& map,
                   int& counter,
                   int max );
};

} // inciter::

#endif // Conductor_h
