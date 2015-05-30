//******************************************************************************
/*!
  \file      src/Inciter/InciterSetup.C
  \author    J. Bakosi
  \date      Sat 30 May 2015 10:46:27 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Functions used to setup inciter
  \details   Functions used to setup inciter.
*/
//******************************************************************************

#include <set>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <type_traits>
#include <algorithm>

#include <boost/format.hpp>
#include <boost/optional.hpp>

#include "Types.h"
#include "Print.h"
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
#include "Inciter/CmdLine/Parser.h"
#include "Inciter/InputDeck/Parser.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "mpi.h"
#include "mpi-interoperate.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

//! Error message to output when the overdecomposition is too large
std::string over = R"(Overdecomposition of the mesh is too large compared to the
number of work units computed based on the degree of virtualization desired. As
a result, there would be at least one work unit with no mesh elements to work
on, i.e., nothing to do. Solution 1: decrease the virtualization to a lower
value using the command-line argument '-u'. Solution 2: decrease the number
processing elements (PEs) using the charmrun command-line argument '+pN' where N
is the number of PEs, which implicitlyincreases the size (and thus decreases the
number) of work units.)";

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
meshinfo( const tk::Print& print,
          const tk::UnsMesh& graph,
          uint64_t load,
          uint64_t chunksize,
          uint64_t remainder,
          uint64_t nchare,
          tk::real virtualization )
//******************************************************************************
//! Print information on mesh graph and load distribution
//! \param[in] print Pretty printer
//! \param[in] graph Unstructured mesh graph object
//! \param[in] load Computational load (here: number of graph nodes)
//! \param[in] chunksize Chunk size, see Base/LoadDistribution.h
//! \param[in] remainder Remainder, see Base/LoadDistribution.h
//! \param[in] nchare Number of work units (Charm++ chares)
//! \param[in] virtualization Degree of virtualization
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
    print.item( "Virtualization [0.0...1.0]", virtualization );
    print.item( "Load (number of tetrahedra)", load );
    print.item( "Number of processing elements", numpes );
    print.item( "Number of work units",
                std::to_string( nchare ) + " (" +
                std::to_string( nchare-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );
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

static std::vector< std::vector< std::size_t > >
elemOwner( std::size_t nchare, const std::vector< std::size_t >& che )
//******************************************************************************
//! Construct global mesh element ids for each chare
//! \param[in] nchare Number chares points are partitioned into
//! \param[in] che Chares of elements: array of chare ownership IDs mapping
//!   graph elements to Charm++ chares. Size: number of elements in mesh graph.
//! \return Vectors of global mesh element ids owned by each chare
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< std::vector< std::size_t > > element( nchare );

  for (std::size_t i=0; i<nchare; ++i)              // for all colors
    for (std::size_t e=0; e<che.size(); ++e )       // for all mesh elements
      if (che[e] == i)
        element[i].push_back( e );

  // This check should always be done, as it can result from incorrect user
  // input compared to the mesh size and not due to programmer error.
  for(const auto& c : element) ErrChk( !c.empty(), std::move(over) );

  return element;
}

static std::vector< std::size_t >
chElem( std::size_t nchare,
        const std::vector< std::size_t >& chp,
        const std::vector< std::size_t >& tetinpoel,
        const std::pair< std::vector< std::size_t >,
                         std::vector< std::size_t > >& esup )
//******************************************************************************
//! Construct array of chare ownership IDs mapping mesh elements to chares
//! \param[in] nchare Number chares points are partitioned into
//! \param[in] chp Chares of points: array of chare ownership IDs mapping graph
//!   points to Charm++ chares. Size: number of points in mesh graph.
//! \param[inout] tetinpoel Mesh tetrahedron element connectivity
//! \param[inout] esup Derived data structure, storing elements surrounding
//!   points in mesh
//! \details This function constructs the vector che, 'chares of elements'. This
//!   is the equivalent of chp, 'chares of points', but stores the chare ids of
//!   mesh elements. The element ownership is computed based on the point
//!   ownership. If all points of an element are owned by a chare, then that
//!   chare owns the element. If not all points of an element are owned by the
//!   a chare, than the chare with lower chare id gets to own the element.
//! \author J. Bakosi
//******************************************************************************
{
  // Lambda to find out if all 4 points of tetrahedron e are owned by the same
  // chare that owns point p; if so, return true, if not, return the lowest of
  // the owner chare ids of the points of the tetrahedron (chosen as the owner
  // of e)
  auto own = [ &chp, &tetinpoel ]( std::size_t e, std::size_t p )
           -> std::pair< bool, std::size_t >
  {
    std::vector< bool > op;
    std::set< std::size_t > owner;
    for (std::size_t n=0; n<4; ++n)
      if (chp[ tetinpoel[e*4+n] ] == chp[p])
        op.push_back( true );
      else
        owner.insert( chp[ tetinpoel[e*4+n] ] );
    if (op.size() == 4)
      return { true, 0 };
    else
      return { false, std::min( chp[p], *owner.begin() ) };
  };

  // Construct array of chare ownership IDs mapping mesh elements to chares
  std::vector< std::size_t > che( tetinpoel.size()/4 );
  for (std::size_t p=0; p<chp.size(); ++p)  // for all mesh points
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

  Assert( nchare > *std::max_element( begin(che), end(che) ),
          "Elements assigned to more than the number chares" );

  // This check should always be done, as it can result from incorrect user
  // input compared to the mesh size and not due to programmer error.
  auto minmax = std::minmax_element( begin(che), end(che) );
  ErrChk( *minmax.first == 0 && *minmax.second == nchare-1, std::move(over) );

  return che;
}

static void
assignMesh(
  const tk::UnsMesh& graph,
  const std::pair< std::vector< std::size_t >,
                   std::vector< std::size_t > >& part,
  const ctr::InputDeck& inputdeck,
  std::size_t& npoin,
  std::vector< std::vector< std::size_t > >& point,
  std::vector< std::vector< std::size_t > >& element,
  std::vector< std::size_t >& meshfilemap,
  std::vector< std::size_t >& tetinpoel,
  std::pair< std::vector< std::size_t >, std::vector< std::size_t > >& esup,
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >& pcomm,
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >& ecomm )
//******************************************************************************
//! Assign mesh partitions to Charm++ chares
//! \param[in] graph Unstructured mesh graph object reference
//! \param[in] part Array of chare ownership IDs mapping graph points to
//!   Charm++ chares, and new->old mesh point id map (new: renumbered, old: as
//!   in mesh file).
//! \param[inout] inputdeck Input deck object filled during parsing user input
//! \param[inout] npoin Total number of points in mesh
//! \param[inout] point Global mesh point ids owned by each chare
//! \param[inout] element Global mesh element ids owned by each chare
//! \param[inout] meshfilemap Index map between renumbered mesh point ids and
//!   those in mesh file
//! \param[inout] tetinpoel Mesh tetrahedron element connectivity
//! \param[inout] esup Derived data structure, storing elements surrounding
//!   points in mesh
//! \param[inout] pcomm Vector of point-based export maps for all chares (empty
//!   if nchares=1)
//! \param[inout] ecomm Vector of elem-based export maps for all chares (empty
//!   if nchares=1)
//! \note This function only operates on MPI rank 0, since that is the only rank
//!   where the input argument 'part' is expected to be non-empty.
//! \author J. Bakosi
//******************************************************************************
{
  int peid;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );

  if (peid == 0) {

    const auto& chp = part.first;

    Assert( chp.size() == graph.size(),
            "Size of ownership array must equal the number of graph nodes" );

    // Return total number points
    npoin = chp.size();

    // Find out number of chares desired
    auto minmax = std::minmax_element( begin(chp), end(chp) );
    auto nchare = *minmax.second - *minmax.first + 1;

    // Construct and return global mesh point ids for each chare
    point = poinOwner( nchare, chp );

    // Return index map between renumbered mesh points and those in mesh file
    meshfilemap = part.second;

    // Return element connectivity
    tetinpoel = graph.tetinpoel();

    // Generate and return elements surrounding points
    esup = tk::genEsup( tetinpoel, 4 );

    //! Construct array of chare ownership IDs mapping mesh elements to chares
    const auto che = chElem( nchare, chp, tetinpoel, esup );

    // Construct and return global mesh point ids for each chare
    element = elemOwner( nchare, che );

    // If graph not partitioned, quit leaving communication maps empty
    if (nchare == 1) return;

    // Compute and return point-based communication maps
    pcomm = tk::poinCommMaps( graph, chp, tetinpoel, nchare, std::move(over) );

    // Compute and return element-based communication maps
    ecomm = tk::elemCommMaps( chp, tetinpoel, element, nchare );
 }
}

ctr::CmdLine
parseCmdLine( int argc, char** argv )
//******************************************************************************
//  Parse command line
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
  ctr::CmdLine cmdline;

  // Parse command line into cmdline
  CmdLineParser cmdParser( argc, argv, print, cmdline, peid );

  // Echo help if requested and make sure mandatory command-line args are set
  help( argc, argv, cmdline, print );

  return cmdline;
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

  InputDeckParser inputdeckParser( print, cmdline, inputdeck, peid );

  if (peid == 0) {
    print.item( "Parsed control file", "success" );
    print.endpart();
    print.part( "Factory" );
  }
}

void
prepareMesh(
  const ctr::CmdLine& cmdline,
  const tk::Print& print,
  const ctr::InputDeck& inputdeck,
  std::vector< std::pair< std::string, tk::Timer::Watch > >& timestamp,
  std::size_t& npoin,
  std::vector< std::vector< std::size_t > >& point,
  std::vector< std::vector< std::size_t > >& element,
  std::vector< std::size_t >& meshfilemap,
  std::vector< std::size_t >& tetinpoel,
  std::pair< std::vector< std::size_t >, std::vector< std::size_t > >& esup,
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >& pcomm,
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >& ecomm )
//******************************************************************************
//! Prepare computational mesh
//! \param[in] cmdline Command-line stack
//! \param[in] print Pretty printer
//! \param[in] inputdeck Input deck with inser input
//! \param[inout] timestamp Time stamps in h:m:s format
//! \param[inout] npoin Total number of points in mesh
//! \param[inout] point Global mesh point ids owned by each chare
//! \param[inout] element Global mesh element ids owned by each chare
//! \param[inout] meshfilemap Index map between renumbered mesh point ids and
//!   those in mesh file
//! \param[inout] tetinpoel Mesh tetrahedron element connectivity
//! \param[inout] esup Derived data structure, storing elements surrounding
//!   points in mesh
//! \param[inout] pcomm Vector of point-based export maps for all chares (empty
//!   if nchares=1)
//! \param[inout] ecomm Vector of elem-based export maps for all chares (empty
//!   if nchares=1)
//! \author J. Bakosi
//******************************************************************************
{
  // The load is taken to be proportional to the number of elements of the
  // mesh which is proportional to the number of points as well as the number
  // of unique edges in the mesh. For a typical mesh of tetrahedra
  //  nelem = 5.5*npoin
  //  nedge = 7*npoin
  //  npsup = 14*npoin,
  // where
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

  int peid, numpes;
  MPI_Status status;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  // Read mesh graph from file only on MPI rank 0 and distribute load size
  if (peid == 0) {
    tk::Timer timer;
    tk::ExodusIIMeshReader er( cmdline.get< tag::io, tag::input >(), graph );
    er.readGraph();
    timestamp.emplace_back( "Read mesh graph from file", timer.hms() );
    load = graph.tetinpoel().size()/4;
    for (int i=1; i<numpes; ++i)
      MPI_Send( &load, 1, MPI_UINT64_T, i, load_tag, MPI_COMM_WORLD );
  } else {
    MPI_Recv( &load, 1, MPI_UINT64_T, 0, load_tag, MPI_COMM_WORLD, &status );
  }

  // Compute load distribution given total work (load) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  const auto nchare =
    tk::linearLoadDistributor( cmdline.get< tag::virtualization >(),
                               load, numpes, chunksize, remainder );

  // Print out info on mesh and load distribution
  meshinfo( print, graph, load, chunksize, remainder, nchare,
            cmdline.get< tag::virtualization >() );

  // Partition graph using Zoltan and assign mesh partitions to Charm++ chares
  tk::Timer t;
  assignMesh( graph, tk::zoltan::partitionMesh(graph, nchare, print),
              inputdeck, npoin, point, element, meshfilemap, tetinpoel, esup,
              pcomm, ecomm );
  timestamp.emplace_back( "Partition mesh & assign to Charm++ chares",
                           t.hms() );
}

} // inciter::
