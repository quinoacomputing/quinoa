//******************************************************************************
/*!
  \file      src/Inciter/InciterSetup.C
  \author    J. Bakosi
  \date      Sat 21 Nov 2015 02:59:23 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Functions used to setup inciter
  \details   Functions used to setup inciter.
*/
//******************************************************************************

#include <string>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <unordered_map>

#include <boost/format.hpp>
#include <boost/optional.hpp>

#include "Config.h"
#include "Types.h"
#include "InciterPrint.h"
#include "Init.h"
#include "Tags.h"
#include "Exception.h"
#include "ExceptionMPI.h"
#include "Keywords.h"
#include "HelpFactory.h"
#include "UnsMesh.h"
#include "InciterSetup.h"
#include "ExodusIIMeshReader.h"
#include "LoadDistributor.h"
#include "ZoltanInterOp.h"
#include "DerivedData.h"
#include "CommMap.h"
#include "Reorder.h"
#include "Inciter/CmdLine/CmdLine.h"
#include "Inciter/CmdLine/Parser.h"
#include "Inciter/InputDeck/Parser.h"
#include "Options/PartitioningAlgorithm.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <mpi.h>
#include <mpi-interoperate.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

//! \brief Error message to output when the overdecomposition is too large
//!   compared to the number of total cells in the mesh
std::string over = "Overdecomposition of the mesh is too large compared to the "
"number of work units computed based on the degree of virtualization desired. "
"As a result, there would be at least one work unit with no mesh elements to "
"work on, i.e., nothing to do. Solution 1: decrease the virtualization to a "
"lower value using the command-line argument '-u'. Solution 2: decrease the "
"number processing elements (PEs) using the charmrun command-line argument "
"'+pN' where N is the number of PEs, which implicitlyincreases the size (and "
"thus decreases the number) of work units.";

static void
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

static void
meshinfo( const ctr::InputDeck& inputdeck,
          const InciterPrint& print,
          const tk::UnsMesh& graph,
          std::size_t npoin,
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
//! \return Total number of elements in mesh (aggregated across all MPI ranks)
//! \author J. Bakosi
//******************************************************************************
{
  int peid, numpes;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  std::size_t penelem = graph.nelem();
  std::size_t nelem = 0;
  MPI_Reduce( &penelem, &nelem, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

  std::size_t pels = graph.lininpoel().size()/2;
  std::size_t ls = 0;
  MPI_Reduce( &pels, &ls, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

  std::size_t pers = graph.triinpoel().size()/3;
  std::size_t rs = 0;
  MPI_Reduce( &pers, &rs, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

  std::size_t pets = graph.tetinpoel().size()/4;
  std::size_t ts = 0;
  MPI_Reduce( &pets, &ts, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

  if (peid == 0) {
    // Print out mesh graph stats
    print.section( "Input mesh graph statistics" );
    print.item( "Number of element blocks", graph.neblk() );
    print.item( "Number of elements", nelem );
    print.item( "Number of nodes", npoin );

    if (ls) print.item( "Number of lines", ls );
    if (rs) print.item( "Number of triangles", rs );
    if (ts) print.item( "Number of tetrahedra", ts );

    // Print out info on load distribution
    print.section( "Load distribution" );
    print.item( "Virtualization [0.0...1.0]",
                inputdeck.get< tag::cmd, tag::virtualization >() );
    print.item( "Load (number of tetrahedra)", load );
    print.item( "Number of processing elements", numpes );
    print.item( "Number of work units",
                std::to_string( nchare ) + " (" +
                std::to_string( nchare-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );
    print.endsubsection();

    // Print out mesh partitioning configuration
    print.section( "Initial mesh partitioning" );
    print.Item< tk::ctr::PartitioningAlgorithm,
                tag::selected, tag::partitioner >();
    print.endsubsection();
  }
}

static std::vector< std::vector< std::size_t > >
poinOwner( std::size_t nchare, const std::vector< std::size_t >& chp )
//******************************************************************************
//! Construct global mesh point ids for each chare
//! \param[in] nchare Number chares points are partitioned into
//! \param[in] chp Chares of points: array of chare ownership IDs mapping graph
//!   points to Charm++ chares. Size: number of points in mesh graph.
//! \return Vectors of global mesh point ids owned by each chare
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< std::vector< std::size_t > > point( nchare );

  for (std::size_t i=0; i<nchare; ++i)              // for all colors
    for (std::size_t p=0; p<chp.size(); ++p )       // for all mesh points
      if (chp[p] == i)
        point[i].push_back( p );

  // This check should always be done, as it can result from incorrect user
  // input compared to the mesh size and not due to programmer error.
  for(const auto& c : point) ErrChk( !c.empty(), std::move(over) );

  return point;
}

static std::unordered_map< int, std::vector< std::size_t > >
elemOwner( const std::vector< std::size_t >& che,
           const std::vector< std::size_t >& geid )
//******************************************************************************
//! Construct global mesh element ids for each chare
//! \param[in] che Chares of elements: array of chare ownership IDs mapping
//!   graph elements to Charm++ chares. Size: number of elements in mesh graph.
//! \param[in] geid List of global element ids for this MPI rank
//! \return Vector of global mesh element ids owned by each chare on this rank
//! \note This function operates on a chunk of the total mesh, distributed among
//!   all MPI ranks.
//! \author J. Bakosi
//******************************************************************************
{
  Assert( che.size() == geid.size(), "The size of the global element index and "
          "the chare element arrays must equal" );

  std::unordered_map< int, std::vector< std::size_t > > element;

  for (std::size_t e=0; e<che.size(); ++e)
    element[ static_cast<int>(che[e]) ].push_back( geid[e] );

  Assert( !element.empty(),
          "No elements assigned to chares on one of the MPI ranks" );

  // This check should always be done, as it can result from incorrect user
  // input compared to the mesh size and not due to programmer error.
  for(const auto& c : element) ErrChk( !c.second.empty(), std::move(over) );

  return element;
}

static std::vector< std::size_t >
chElem( const std::vector< std::size_t >& chp,
        const std::vector< std::size_t >& tetinpoel )
//******************************************************************************
//! Construct array of chare ownership IDs mapping mesh elements to chares
//! \param[in] chp Chares of points: array of chare ownership IDs mapping graph
//!   points to Charm++ chares. Size: number of points in chunk of mesh graph.
//! \param[in] tetinpoel Mesh tetrahedron element connectivity with global IDs
//! \return Vector of chare IDs mapped to mesh elements
//! \details This function constructs the vector che, 'chares of elements'. This
//!   is the equivalent of chp, 'chares of points', but stores the chare ids of
//!   mesh elements. The element ownership is computed based on the point
//!   ownership. If all points of an element are owned by a chare, then that
//!   chare owns the element. If not all points of an element are owned by the
//!   a chare, than the chare with lower chare id gets to own the element.
//! \note This function operates on all MPI ranks, working on different chunks
//!   of the mesh.
//! \author J. Bakosi
//******************************************************************************
{
  // Generate element connectivity storing local node ids
  std::vector< std::size_t > inpoel, gid;
  std::tie( inpoel, gid ) = tk::global2local( tetinpoel );

  // Lambda to find out if all 4 points of tetrahedron e are owned by the same
  // chare that owns point p; if so, return true, if not, return the lowest of
  // the owner chare ids of the points of the tetrahedron (chosen as the owner
  // of e)
  auto own = [ &chp, &inpoel ]( std::size_t e, std::size_t p )
           -> std::pair< bool, std::size_t >
  {
    std::vector< bool > op;
    std::set< std::size_t > owner;
    for (std::size_t n=0; n<4; ++n)
      if (chp[ inpoel[e*4+n] ] == chp[p])
        op.push_back( true );
      else
        owner.insert( chp[ inpoel[e*4+n] ] );
    if (op.size() == 4)
      return { true, 0 };
    else
      return { false, std::min( chp[p], *owner.begin() ) };
  };

  // Generate elements surrounding points based on connectivity
  auto esup = tk::genEsup( inpoel, 4 );

  // Construct array of chare ownership IDs mapping mesh elements to chares
  std::vector< std::size_t > che( inpoel.size()/4 );
  for (std::size_t p=0; p<chp.size(); ++p)  // for all mesh points in this chunk
    for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i) {
      auto e = esup.first[i];
      auto o = own(e,p);
      if (o.first)
        che[e] = chp[p];
      else
        che[e] = o.second;
    }

//   std::cout << "che: ";
//   for (auto o : che) std::cout << o << " ";
//   std::cout << '\n';
// 
//   std::cout << "chp: ";
//   for (auto o : chp) std::cout << o << " ";
//   std::cout << '\n';

  return che;
}

void
parseCmdLine( int argc, char** argv, ctr::CmdLine& cmdline )
//******************************************************************************
//  Parse command line
//! \param[in] argc Number of command-line arguments passed to executable
//! \param[in] argv C-style character array of character arrays (cmd line args)
//! \param[inout] cmdline Command line data to fill by command line parser
//! \return Hierarchical tagged tuple containind command line data
//! \author J. Bakosi
//******************************************************************************
{
  int peid;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );

  // Create basic pretty printer
  tk::Print print;    // quiet output by default using print, see tk::Print ctor

  // Parse command line into cmdline
  CmdLineParser cmdParser( argc, argv, print, cmdline, peid );

  // Echo help if requested and make sure mandatory command-line args are set
  help( argc, argv, cmdline, print );
}

void
init( const ctr::CmdLine& cmdline,
      const tk::Print& print,
      ctr::InputDeck& inputdeck,
      int argc,
      char** argv )
//******************************************************************************
//! Parse command line and input deck, instantiate pretty printer, echo info
//! \param[in] cmdline Command-line stack
//! \param[in] print Pretty printer
//! \param[inout] inputdeck Input deck object filled during parsing user input
//! \param[in] argc Number of command-line arguments passed to executable
//! \param[in] argv C-style character array of character arrays (cmd line args)
//! \author J. Bakosi
//******************************************************************************
{
  int peid;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );

  if (peid == 0) {
    // Echo program header
    echoHeader( print, tk::HeaderType::INCITER );
    // Echo environment
    print.part( "Environment" );
    // Build environment
    echoBuildEnv( print, INCITER_EXECUTABLE );
    // Runtime environment
    echoRunEnv( print, argc, argv, cmdline.get< tag::verbose >() );
  }

  // Parse input deck into inputdeck
  if (peid == 0)
    print.item( "Control file", cmdline.get< tag::io, tag::control >() );

  // Parse input deck, also store parsed command line in input deck
  InputDeckParser( print, cmdline, inputdeck, peid );

  if (peid == 0) {
    print.item( "Parsed control file", "success" );
    print.endpart();
    print.part( "Factory" );
  }
}

std::unordered_map< int, std::vector< std::size_t > >
prepareMesh(
  const ctr::InputDeck& inputdeck,
  const InciterPrint& print,
  std::vector< std::pair< std::string, tk::Timer::Watch > >& timestamp,
  uint64_t& nchare,
  std::size_t& npoin )
//******************************************************************************
//! Prepare computational mesh
//! \param[in] inputdeck Input deck with inser input
//! \param[in] print Pretty printer
//! \param[inout] timestamp Time stamps in h:m:s format
//! \param[inout] npoin Total number of points in mesh
//! \return Global mesh element ids owned by each chare
//! \details The load is taken to be proportional to the number of elements of
//!   the mesh which is proportional to the number of points as well as the
//!   number of unique edges in the mesh. For a typical mesh of tetrahedra
//!   - nelem = 5.5*npoin
//!   - nedge = 7*npoin
//!   - npsup = 14*npoin,
//!   where
//!   - nelem - number of elements,
//!   - npoin - number of points,
//!   - nedge - number of unique edges,
//!   - npsup - number of points surrounding points, which is the same as the
//!  number of (non-unique) edges surrounding points.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008.
//! \author J. Bakosi
//******************************************************************************
{
  int peid, numpes;
  MPI_Status status;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  // Read mesh graph from file, a chunk on each MPI rank 
  tk::Timer timer;
  if (peid == 0)
    print.diagstart( "Reading mesh graph ..." );

  const auto& cmdline = inputdeck.get< tag::cmd >();

  // Get number of mesh points and number of tetrahedron elements in input file
  tk::ExodusIIMeshReader er( cmdline.get< tag::io, tag::input >() );
  npoin = er.readElemBlockIDs();
  auto nel = er.nel( tk::ExoElemType::TET );

  // Read our chunk of tetrahedron element connectivity from file and also store
  // the list of global element indices for our chunk of the mesh
  auto chunk = nel / numpes;
  auto from = peid * chunk;
  auto till = from+chunk;
  if (peid == numpes-1) till += nel % numpes;
  std::vector< std::size_t > tetinpoel;
  std::array< std::size_t, 2 > ext = { { static_cast<std::size_t>(from),
                                         static_cast<std::size_t>(till) } };
  er.readElements( ext, tk::ExoElemType::TET, tetinpoel );
  std::vector< std::size_t > geid( static_cast<std::size_t>(till-from) );
  std::iota( begin(geid), end(geid), from );

  if (peid == 0) print.diagend( "done" );

  // If geometric partitioner is selected, compute element centroid coordinates
  std::array< std::vector< tk::real >, 3 > centroid;
  const auto algorithm = inputdeck.get< tag::selected, tag::partitioner >();
  if ( algorithm == tk::ctr::PartitioningAlgorithmType::RCB ||
       algorithm == tk::ctr::PartitioningAlgorithmType::RIB )
  {
    if (peid == 0)
      print.diagstart("Reading mesh nodes and computing element centroids ...");
    // Collect global node ids of elements of this rank
    auto gid = tetinpoel;
    tk::unique( gid );
    // Read node coordinates of mesh elements (of this rank) from file
    auto coord = er.readNodes( tk::extents(gid) );
    // Make room for element centroid coordinates
    auto& cx = centroid[0];
    auto& cy = centroid[1];
    auto& cz = centroid[2];
    auto num = tetinpoel.size()/4;
    cx.resize( num );
    cy.resize( num );
    cz.resize( num );
    // Compute element centroids for our chunk of the mesh elements
    for (std::size_t e=0; e<num; ++e) {
      const auto& a = tk::cref_find( coord, tetinpoel[e*4+0] );
      const auto& b = tk::cref_find( coord, tetinpoel[e*4+1] );
      const auto& c = tk::cref_find( coord, tetinpoel[e*4+2] );
      const auto& d = tk::cref_find( coord, tetinpoel[e*4+3] );
      cx[e] = (a[0] + b[0] + c[0] + d[0]) / 4.0;
      cy[e] = (a[1] + b[1] + c[1] + d[1]) / 4.0;
      cz[e] = (a[2] + b[2] + c[2] + d[2]) / 4.0;
    }
    if (peid == 0) print.diagend( "done" );
  }

  // Create empty unstructured mesh object initializing only connectivity, so we
  // call it a graph instead of a mesh; note that tetinpoel is moved
  tk::UnsMesh graph( std::move( tetinpoel ) );

  // Compute and distribute total load to all PEs
  uint64_t peload = graph.tetinpoel().size()/4;
  uint64_t load = 0;
  MPI_Allreduce( &peload, &load, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

  timestamp.emplace_back( "Read mesh from file", timer.hms() );

  // Compute load distribution given total work (load) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  nchare = tk::linearLoadDistributor( cmdline.get< tag::virtualization >(),
                                      load, numpes, chunksize, remainder );

  // Print out info on mesh and load distribution
  meshinfo( inputdeck, print, graph, npoin, load, chunksize, remainder, nchare );

  // Partition graph using Zoltan and assign mesh partitions to Charm++ chares
  tk::Timer t;
  if (peid == 0)
    print.diagstart( "Partitioning mesh & assigning to Charm++ chares..." );

  const auto che =
    tk::zoltan::geomPartMesh( algorithm, centroid, geid,
                              graph.tetinpoel().size()/4, nchare );

  Assert( che.size() == graph.tetinpoel().size()/4, "Size of ownership array "
          "does not equal the number of graph elements" );

  // Construct and return global mesh element ids for each chare
  auto element = elemOwner( che, geid );

  if (peid == 0) print.diagend( "done" );
  timestamp.emplace_back( "Partition mesh & assign to Charm++ chares",
                           t.hms() );

  return element;
}

} // inciter::
