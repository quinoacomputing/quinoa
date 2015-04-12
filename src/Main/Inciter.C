//******************************************************************************
/*!
  \file      src/Main/Inciter.C
  \author    J. Bakosi
  \date      Sun 12 Apr 2015 07:22:27 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter, computational shock hydrodynamics tool, Charm++ main
    chare.
  \details   Inciter, computational shock hydrodynamics tool, Charm++ main
    chare. This file contains the definition of the Charm++ main chare,
    equivalent to main() in Charm++-land.
*/
//******************************************************************************
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <mpi.h>
#include <mpi-interoperate.h>   // for interoperation of MPI and Charm++

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <pup_stl.h>

#include <Config.h>
#include <RNG.h>
#include <RNGStack.h>
#include <InciterPrint.h>
#include <InciterDriver.h>
#include <Inciter/CmdLine/Parser.h>
#include <Inciter/InputDeck/Parser.h>
#include <MeshFactory.h>
#include <LoadDistributor.h>
#include <ZoltanInterOp.h>
#include <ExceptionMPI.h>
#include <ProcessException.h>
#include <ExodusIIMeshReader.h>
#include <inciter.decl.h>
#include <Init.h>

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

//! Inciter declarations and definitions
namespace inciter {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs
//! by the Charm++ runtime system. Though semantically not const, all these
//! global data should be considered read-only. See also
//! http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. The data
//! below is global-scope because they must be available to all PEs which could
//! be on different machines.

//! Defaults of input deck, facilitates detection what is set by user
ctr::InputDeck g_inputdeck_defaults;
//! Input deck filled by parser, containing all input data
ctr::InputDeck g_inputdeck;
//! Derived data structure, storing elements surrounding points in mesh
std::pair< std::vector< std::size_t >, std::vector< std::size_t > > g_esup;
//! Mesh tetrahedron element connectivity
std::vector< int > g_tetinpoel;
//! Graph coloring for all mesh points
std::vector< std::size_t > g_colors;
//! Vector of export/import maps for all chare ids (empty if nchares = 1)
std::vector< std::map< std::size_t, std::vector< std::size_t > > > g_comm;

//! \brief Time stamps in h:m:s for the initial MPI portion
//! \details Time stamps collected here are those collected by the initial MPI
//!   portion and are displayed by the Charm++ main chare at the end. While this
//!   map of timers is declared in global scope (so that the Charm++ main chare
//!   can access it), it is intentionally NOT declared in the Charm++ main
//!   module interface file for Inciter in Main/inciter.ci, so that the Charm++
//!   runtime system does not migrate it across all PEs. This is okay, since
//!   since there is no need for any of the other Charm++ chares to access it in
//!   the future. In fact, the main chare grabs it and swallows it right away
//!   during its constructor.
std::vector< std::pair< std::string, tk::Timer::Watch > > g_timestamp;

//! Conductor Charm++ proxy facilitating call-back to Conductor by the
//! individual performers
CProxy_Conductor g_ConductorProxy;

} // inciter::

//! \brief Charm++ main chare for the shock hydroddynamics executable, inciter.
//! \details In inciter the Charm++ runtime system is initialized only after the
//!   mesh has been read in, partitioned, and the necessary data structures,
//!   e.g., communication maps, have been generated. This delayed initialization
//!   of the Charm++ runtime system is required since the mesh partitioning is
//!   done by Zoltan, an MPI library. Note that this Charm++ main chare object
//!   should not be in a namespace.
//! \author J. Bakosi
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details The main chare constructor is the main entry point of the
    //!   Charm++ portion of inciter, called by the Charm++ runtime system. The
    //!   constructor does basic initialization steps, prints out some useful
    //!   information to screen (in verbose mode), and instantiates a driver.
    //!   Since Charm++ is fully asynchronous, the constructor usually spawns
    //!   asynchronous objects and immediately exits. Thus in the body of the
    //!   main chare constructor we fire up an 'execute' chare, which then calls
    //!   back to Main::execute(). Finishing the main chare constructor the
    //!   Charm++ runtime system then starts the network-migration of all
    //!   global-scope data (if any). The execute chare calling back to
    //!   Main::execute() signals the end of the migration of the global-scope
    //!   data. Then we are ready to execute the driver. Since inciter is
    //!   parallel and asynchronous, its driver fires up additional Charm++
    //!   chare objects which then call back to Main::finalize() at some point
    //!   in the future when all work has been finished. finalize() then exits
    //!   by calling Charm++'s CkExit(), shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    Main( CkArgMsg* msg )
    try :
      // Create pretty printer initializing output streams based on command line
      m_print( inciter::g_inputdeck.get< tag::cmd, tag::verbose >()
                 ? std::cout : std::clog ),
      // Create Inciter driver
      m_driver( inciter::InciterDriver( m_print ) ),
      // Start new timer measuring the total runtime
      m_timer(1),
      // Import, i.e., swallow, timers from the initial MPI portion
      m_timestamp( std::move(inciter::g_timestamp) )
    {
      const auto& cmdline = inciter::g_inputdeck.get< tag::cmd >();
      const auto helpcmd = cmdline.get< tag::help >();
      const auto helpctr = cmdline.get< tag::helpctr >();
      const auto helpkw = cmdline.get< tag::helpkw >();
      // Exit if help was requested or exectuable was called without argument
      if (msg->argc == 1 || helpcmd || helpctr || !helpkw.keyword.empty()) {
        CkExit();
      } else {  // business as usual
        mainProxy = thisProxy;
        // Fire up an asynchronous execute object, which when created at some
        // future point in time will call back to this->execute(). This is
        // necessary so that this->execute() can access already migrated
        // global-scope data.
        CProxy_execute::ckNew();
        // Start new timer measuring the migration of global-scope data
        m_timer.emplace_back();
      }
      delete msg;
    } catch (...) { tk::processExceptionCharm(); }

    //! Execute driver created and initialized by constructor
    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        m_driver.execute();
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Normal exit point
    void finalize() {
      try {
        if (!m_timer.empty()) {
          m_timestamp.emplace_back("Total Charm++ runtime", m_timer[0].hms());
          m_print.time( "Timers (h:m:s)", m_timestamp );
          m_print.endpart();
        }
      } catch (...) { tk::processExceptionCharm(); }
      // Tell the Charm++ runtime system to exit
      CkExit();
    }

    //! Add time stamp contributing to final timers output
    void timestamp( std::string label, tk::real stamp ) {
      try {
        m_timestamp.emplace_back( label, tk::hms( stamp ) );
      } catch (...) { tk::processExceptionCharm(); }
    }

  private:
    inciter::InciterPrint m_print;                    //!< Pretty printer
    inciter::InciterDriver m_driver;                  //!< Driver
    std::vector< tk::Timer > m_timer;                 //!< Timers

    //! Time stamps in h:m:s with labels
    std::vector< std::pair< std::string, tk::Timer::Watch > > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
//! \author J. Bakosi
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

namespace inciter {

void
help( int argc,
      char** argv,
      const ctr::CmdLine& cmdline,
      const tk::Print& print )
//******************************************************************************
//! Echo help if requested and make sure mandatory command-line args are set
//! \param[in] argc Number of command-line arguments passed to executable
//! \param[in] argv C-style character array of character arrays (cmd line args)
//! \param[in] cmdline Command-line stack
//! \param[in] print Pretty printer
//! \author J. Bakosi
//******************************************************************************
{
  int peid;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );

  const auto helpcmd = cmdline.get< tag::help >();
  const auto helpctr = cmdline.get< tag::helpctr >();
  const auto helpkw = cmdline.get< tag::helpkw >();

  // Output help if requested or exectuable was called without argument and exit
  if (argc == 1 || helpcmd || helpctr || !helpkw.keyword.empty()) {

    // Initialize the Charm++ runtime system to output Charm++/Converse help
    CharmLibInit( MPI_COMM_WORLD, argc, argv );
    CharmLibExit();

    // Print out help on all command-line arguments if the executable was
    // invoked without arguments or help was requested
    if (peid == 0 && (argc == 1 || helpcmd))
      print.help< tk::QUIET >( INCITER_EXECUTABLE,
                               cmdline.get< tag::cmdinfo >(),
                               "Command-line Parameters:", "-" );

    // Print out help on all control file keywords if they were requested
    if (peid == 0 && helpctr)
      print.help< tk::QUIET >( INCITER_EXECUTABLE,
                               cmdline.get< tag::ctrinfo >(),
                               "Control File Keywords:" );

    // Print out verbose help for a single keyword if requested
    if (peid == 0 && !helpkw.keyword.empty())
      print.helpkw< tk::QUIET >( INCITER_EXECUTABLE, helpkw );

    // Quit
    MPI_Finalize();
    exit( tk::ErrCode::SUCCESS );
  }

  // Make sure mandatory command-line arguments are set
  auto ctralias = kw::control().alias();
  ErrChkMPI( !(cmdline.get< tag::io, tag::control >().empty()),
             "Mandatory control file not specified. "
             "Use '--" + kw::control().string() + " <filename>'" +
             ( ctralias ? " or '-" + *ctralias + " <filename>'" : "" ) + '.' );

  auto inpalias = kw::input().alias();
  ErrChkMPI( !(cmdline.get< tag::io, tag::input >().empty()),
             "Mandatory input file not specified. "
             "Use '--" + kw::input().string() + " <filename>'" +
             ( inpalias ? " or '-" + *inpalias + " <filename>'" : "" ) + '.' );
}

inciter::ctr::CmdLine
parseCmdLine( int argc, char** argv )
//******************************************************************************
//! Parse command line
//! \param[in] argc Number of command-line arguments passed to executable
//! \param[in] argv C-style character array of character arrays (cmd line args)
//! \return Hierarchical tagged tuple containind command line data
//! \author J. Bakosi
//******************************************************************************
{
  int peid;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );

  // Create basic pretty printer
  tk::Print print;    // quiet output by default using print, see tk::Print ctor

  // Create hearchical tagged tuple to store data from command line
  inciter::ctr::CmdLine cmdline;

  // Parse command line into cmdline
  inciter::CmdLineParser cmdParser( argc, argv, print, cmdline, peid );

  // Echo help if requested and make sure mandatory command-line args are set
  inciter::help( argc, argv, cmdline, print );

  return cmdline;
}

void
meshinfo( const tk::Print& print,
          const tk::UnsMesh& graph,
          uint64_t load,
          uint64_t chunksize,
          uint64_t remainder,
          uint64_t nchare )
//******************************************************************************
//! Print information on mesh graph and load distribution
//! \param[in] print Pretty printer
//! \param[in] graph Unstructured mesh graph object
//! \param[in] load Computational load (here: number of graph nodes)
//! \param[in] chunksize Chunk size, see Base/LoadDistribution.h
//! \param[in] remainder Remainder, see Base/LoadDistribution.h
//! \param[in] nchare Number of work units (Charm++ chares)
//! \author J. Bakosi
//******************************************************************************
{
  int peid, numpes;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  if (peid == 0) {
    // Print out mesh graph stats
    print.section( "Input mesh graph statistics" );
    print.item( "Number of element blocks", graph.neblk() );
    print.item( "Number of elements", graph.nelem() );
    print.item( "Number of nodes", graph.size() );

    if (!graph.lininpoel().empty())
      print.item( "Number of lines", graph.lininpoel().size()/2 );
    if (!graph.triinpoel().empty())
      print.item( "Number of triangles", graph.triinpoel().size()/3 );
    if (!graph.tetinpoel().empty())
      print.item( "Number of tetrahedra", graph.tetinpoel().size()/4 );

    // Print out info on load distribution
    print.section( "Load distribution" );
    print.item( "Virtualization [0.0...1.0]",
                  g_inputdeck.get< tag::cmd, tag::virtualization >() );
    print.item( "Load (number of mesh points)", load );
    print.item( "Number of processing elements", numpes );
    print.item( "Number of work units",
                std::to_string( nchare ) + " (" +
                std::to_string( nchare-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );
  }
}

void
comMaps( const tk::UnsMesh& graph,
         tk::tuple::tagged_tuple<
           tag::esup,  std::pair< std::vector< std::size_t >,
                                  std::vector< std::size_t > >,
           tag::psup,  std::pair< std::vector< std::size_t >,
                                  std::vector< std::size_t > >,
           tag::owner, std::vector< std::size_t > >&& partitions )
//******************************************************************************
//! Compute communication (export-, and import-) maps for all graph partitions
//! \param[in] graph Unstructured mesh graph object
//! \param[in] partitions Tagged tuple containing elements surrounding points
//!   (at tag::esup), see tk::genEsup(), points surrounding points (at
//!   tag::psup), see tk::genPsup(), and array of chare ownership IDs mapping
//!   graph points to concurrent arsync chares (at tag::owner). Assumed to be
//!   non-emtpy only on MPI rank 0. Note that this tuple is swallowed in and
//!   cannibalized moving some of its parts to global scope.
//! \return Vector of export/import maps associating receiver/sender chare ids
//!   to unique communicated point ids for all chare ids, only on MPI rank 0.
//! \details Compute communication maps: (1) the export map, associating chare
//!   ids to a set of receiver chare ids associated to unique mesh points sent,
//!   and (2) the import map, associating chare ids to a set of sender chare ids
//!   associated to unique mesh points received. In the MPI paradigm, these maps
//!   correspond to the export and import lists, respectively, i.e., lists of
//!   ids exported by a given rank to a set of receiver ranks and their
//!   associated mesh points sent at which data are to be sent (export), and
//!   lists of ids imported by a given rank to a set of sender ranks and their
//!   associated mesh points sent at which data are to be received (import).
//!   For example, in Zoltan, this roughly corresponds to the "exported" and
//!   "imported" ids after partitioning. Actually, Zoltan_LB_Partition() already
//!   returns this information. However, if the partitioning with Zoltan is done
//!   using less MPI ranks than the number of desired mesh partitions, i.e.,
//!   overdecomposition (as is the case here), there are more mesh partitions
//!   than MPI ranks and thus the arrays returned by Zoltan are not sufficient
//!   to determine the export and import maps for all the chares. Thus we
//!   compute the export and import mapping here.
//! \note This function only operates on MPI rank 0, since that is the only rank
//!   where the argument partitions is expected to be non-empty.
//! \author J. Bakosi
//******************************************************************************
{
  int peid;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );

  if (peid == 0) {
    const auto& psup1 = partitions.get< tag::psup >().first;
    const auto& psup2 = partitions.get< tag::psup >().second;
    const auto& owner = partitions.get< tag::owner >();

    Assert( owner.size() == graph.size(),
            "Size of ownership array must equal the number of graph nodes" );

    // Move derived data structure storing elements surrounding points of mesh
    // to global scope for Charm++ chares
    g_esup = std::move( partitions.get< tag::esup >() );

    // find out number of chares desired
    auto minmax = std::minmax_element( begin(owner), end(owner) );
    auto nchare = *minmax.second - *minmax.first + 1;
    // if graph not partitioned, nothing to do, leave communication maps empty
    if (nchare == 1) {
      // Move ownership array to global scope for Charm++ chares
      g_colors = std::move( owner );
      return;
    }

    auto npoin = graph.size();

    // map to associate a chare id to a map of receiver/sender chare ids
    // accociated to unique point ids sent/received (export/import map)
    std::map< std::size_t,
              std::map< std::size_t, std::set< std::size_t > > > comm;

    // construct export and import maps
    for (std::size_t p=0; p<npoin; ++p)
      for (auto i=psup2[p]+1; i<=psup2[p+1]; ++i) {
        auto q = psup1[i];
        if (owner[p] != owner[q])
          comm[ owner[p] ][ owner[q] ].insert( p );
      }

    // This check should always be done, as it can result from incorrect user
    // input compared to the mesh size and not due to programmer error.
    ErrChk( comm.size() == nchare,
            "Number of export/import maps computed (" +
            std::to_string(comm.size()) + ") must equal the number of "
            "work units desired (" + std::to_string(nchare) + "). "
            "This happens when the overdecomposition is too large compared to "
            "the number of work units computed based on the degree of "
            "virtualization desired. As a result, there would be " +
            std::to_string(nchare-comm.size()) +" work unit(s) with nothing to "
            "do. Solution 1: decrease the virtualization (currently: "
            + std::to_string(g_inputdeck.get< tag::cmd, tag::virtualization >())
            + ") to a lower value using the command-line argument '-u'. "
            "Solution 2: decrease the number processing elements (PEs) using "
            "the charmrun command-line argument '+pN' where N is the number of "
            "PEs, which implicitly increases the size (and thus decreases the "
            "number) of work units." );

    std::size_t c = 0;
    for (const auto& e : comm) {
      Assert( e.first == c++,
              "Export/import maps should not be missing for chare id " +
              std::to_string(c-1) );
      for (const auto& x : e.second)
        Assert( x.first >= 0,
                "Export/import map recv/send chare ids must be non-negative" );
    }

    // Move ownership array to global scope for Charm++ chares
    g_colors = std::move( owner );

    // Construct final product: a vector of export/import maps associating
    // receiver/sender chare ids to unique communicated point ids for all chare
    // ids, and store it in global scope so that the main Charm++ chare can
    // access it
    for (const auto& e : comm) {
      g_comm.push_back( {} );
      for (const auto& x : e.second)
        for (auto p : x.second)
          g_comm.back()[ x.first ].push_back( p );
    }

    Assert( g_comm.size() == nchare,
            "Number of export/import maps must equal the number of chares" );

//     std::size_t h = 0;
//     for (const auto& m : g_comm) {
//       std::cout << h++ << " -> ";
//       for (const auto& x : m) {
//         std::cout << x.first << ": ";
//         for (auto p : x.second)
//           std::cout << p << " ";
//       }
//       std::cout << '\n';
//       std::cout << '\n';
//     }
  }
}

} // inciter::

//! \brief Inciter main()
//! \details Inciter does have a main() function so that Zoltan, an MPI library,
//!   can partition the mesh. Then we initialize the Charm++ runtime system and
//!   do a calculation. This is necessary, since MPI_Init() is a bit adamant
//!   about capturing resources it wants and hence it has to be called before
//!   Charm is initialized.
//! \author J. Bakosi
int main( int argc, char **argv ) {

  using inciter::g_timestamp;

  tk::Timer mpi;        // start timing the MPI portion

  int peid, numpes;
  MPI_Status status;

  // Initialize MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  try {

    // Parse command line into cmdline
    const auto cmdline = inciter::parseCmdLine( argc, argv );
   
    // Instantiate inciter's pretty printer
    inciter::InciterPrint
      iprint( cmdline.get< tag::verbose >() ? std::cout : std::clog );

    if (peid == 0) {
      // Echo program header
      echoHeader( iprint, tk::HeaderType::INCITER );
      // Echo environment
      iprint.part( "Environment" );
      // Build environment
      echoBuildEnv( iprint, INCITER_EXECUTABLE );
      // Runtime environment
      echoRunEnv( iprint, argc, argv, cmdline.get< tag::verbose >() );
    }

    // Parse input deck into g_inputdeck
    if (peid == 0)
      iprint.item( "Control file", cmdline.get< tag::io, tag::control >() );

    inciter::InputDeckParser
      inputdeckParser( iprint, cmdline, inciter::g_inputdeck, peid );

    if (peid == 0) {
      iprint.item( "Parsed control file", "success" );
      iprint.endpart();
      iprint.part( "Factory" );
    }

    // The load is taken to be proportional to the number of points of the mesh
    // which is proportional to the number of unique edges in the mesh. Note
    // that for a typical mesh of tetrahedra nelem = 5.5*npoin, nedge = 7*npoin,
    // and npsup = 14*npoin, where
    //  * nelem - number of elements,
    //  * npoin - number of points,
    //  * nedge - number of unique edges,
    //  * npsup - number of points surrounding points, which is the same as the
    // number of (non-unique) edges surrounding points. See also Lohner, An
    // Introduction to Applied CFD Techniques, Wiley, 2008.
    uint64_t load;
    int load_tag = 1;

    // Create empty unstructured mesh object (will only load connectivity, so we
    // call it a graph instead of a mesh)
    tk::UnsMesh graph;

    // Read mesh graph from file only on MPI rank 0 and distribute load size
    if (peid == 0) {
      tk::Timer timer;
      tk::ExodusIIMeshReader er( cmdline.get< tag::io, tag::input >(), graph );
      er.readGraph();
      g_timestamp.emplace_back( "Read mesh graph from file", timer.hms() );
      load = graph.size();
      for (int i=1; i<numpes; ++i)
        MPI_Send( &load, 1, MPI_UINT64_T, i, load_tag, MPI_COMM_WORLD );
    } else {
      MPI_Recv( &load, 1, MPI_UINT64_T, 0, load_tag, MPI_COMM_WORLD, &status );
    }

    // Store tetrahedron element connectivity graph in global scope so Charm++
    // chares will be able to access it
    inciter::g_tetinpoel = graph.tetinpoel();

    // Compute load distribution given total work (load) and user-specified
    // virtualization
    uint64_t chunksize, remainder;
    const auto nchare =
      tk::linearLoadDistributor( cmdline.get< tag::virtualization >(),
                                 load, numpes, chunksize, remainder );

    // Print out info on mesh and load distribution
    inciter::meshinfo( iprint, graph, load, chunksize, remainder, nchare );

    // Partition graph using Zoltan and compute communication maps (stored in
    // g_comMaps) for each graph partition, each of which will become Charm++
    // chares
    tk::Timer t;
    inciter::comMaps( graph, tk::zoltan::partitionMesh(graph, nchare, iprint) );
    g_timestamp.emplace_back("Partition mesh & compute communication maps", t.hms());

  } catch (...) { tk::processExceptionMPI(); }

  g_timestamp.emplace_back( "Total MPI runtime", mpi.hms());

  // Run Charm++ main chare using the partitioned graph
  CharmLibInit( MPI_COMM_WORLD, argc, argv );
  CharmLibExit();

  // Finalize MPI
  MPI_Finalize();

  return tk::ErrCode::SUCCESS;
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <inciter.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
