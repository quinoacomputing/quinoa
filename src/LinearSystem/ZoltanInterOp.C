//******************************************************************************
/*!
  \file      src/Mesh/ZoltanInterOp.C
  \author    J. Bakosi
  \date      Wed 15 Apr 2015 10:07:57 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh graph
    partitioning.
*/
//******************************************************************************
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <zoltan.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <ExceptionMPI.h>
#include <ZoltanInterOp.h>
#include <DerivedData.h>

namespace tk {
namespace zoltan {

//! Zoltan hypergraph data structure
struct HGRAPH_DATA {
  int numMyVertices;            //!< number of vertices that I own initially
  ZOLTAN_ID_TYPE *vtxGID;       //!< global ID of these vertices
  int numMyHEdges;              //!< number of my hyperedges
  int numAllNbors;              //!< number of vertices in my hyperedges
  ZOLTAN_ID_TYPE *edgeGID;      //!< global ID of each of my hyperedges
  int *nborIndex;               //!< index into nborGID array of edge ids
  ZOLTAN_ID_TYPE *nborGID;      //!< array of edge ids
};

static int
get_number_of_vertices( void *data, int *ierr )
//******************************************************************************
//! \brief Zoltan query function returning the number of objects assigned to an
//!   MPI rank
//! \param[inout] data Pointer to user-defined data
//! \param[inout] ierr Error code
//! \return The number of objects that are assigned to the MPI rank
//! \author J. Bakosi
//******************************************************************************
{
  HGRAPH_DATA *hg = (HGRAPH_DATA*)data;
  *ierr = ZOLTAN_OK;
  return hg->numMyVertices;
}

static void get_vertex_list( void *data,
                             int sizeGID,
                             int sizeLID,
                             ZOLTAN_ID_PTR globalID,
                             ZOLTAN_ID_PTR localID,
                             int wgt_dim,
                             float *obj_wgts,
                             int *ierr )
//******************************************************************************
//! \brief Zoltan query function query function to fill two (three if weights
//!    are used) arrays with information about the objects assigned to an MPI
//!    rank
//! \details This Zoltan query function fills two (three if weights are used)
//!   arrays with information about the objects assigned to an MPI rank. Both
//!   arrays are allocated (and subsequently freed) by Zoltan; their size is
//!   determined by a call to a ZOLTAN_NUM_OBJ_FN query function
//!   (get_number_of_vertices) to get the array size. For many algorithms,
//!   either a ZOLTAN_OBJ_LIST_FN query function or a
//!   ZOLTAN_FIRST_OBJ_FN/ZOLTAN_NEXT_OBJ_FN query-function pair must be
//!   registered; however, both query options need not be provided. The
//!   ZOLTAN_OBJ_LIST_FN is preferred for efficiency.
//! \param[inout] data Pointer to user-defined data
//! \param[in] sizeGID The number of array entries used to describe a single
//!   global ID. This value is the maximum value over all MPI ranks of the
//!   parameter NUM_GID_ENTRIES.
//! \param[in] sizeLID The number of array entries used to describe a single
//!   local ID. This value is the maximum value over all MPI ranks of the
//!   parameter NUM_LID_ENTRIES. (It should be zero if local ids are not used.)
//! \param[inout] globalID Upon return, an array of unique global IDs for all
//!   objects assigned to the MPI rank.
//! \param[inout] localID Upon return, an array of local IDs, the meaning of
//!   which can be determined by the application, for all objects assigned to
//!   the MPI rank (optional).
//! \param[in] wgt_dim The number of weights associated with an object
//!   (typically 1), or 0 if weights are not requested. This value is set
//!   through the parameter OBJ_WEIGHT_DIM.
//! \param[inout] obj_wgts Upon return, an array of object weights. Weights for
//!   object i are stored in obj_wgts[(i-1)*wgt_dim:i*wgt_dim-1]. If wgt_dim=0,
//!   the return value of obj_wgts is undefined and may be NULL.
//! \param[inout] ierr Error code
//! \return The number of objects that are assigned to the MPI rank
//! \author J. Bakosi
//******************************************************************************
{
  HGRAPH_DATA *hg = (HGRAPH_DATA*)data;
  *ierr = ZOLTAN_OK;

  // Return the IDs of our vertices, but no weights. Zoltan will assume equally
  // weighted vertices.
  for (int i=0; i<hg->numMyVertices; ++i) {
    globalID[i] = hg->vtxGID[i];
    localID[i] = static_cast< ZOLTAN_ID_TYPE >( i );
  }
}

static void
get_hypergraph_size( void *data,
                     int *num_lists,
                     int *num_nonzeros,
                     int *format,
                     int *ierr )
//******************************************************************************
//! \brief Zoltan query function query function which format we will supply the
//!   hypergraph
//! \details A hypergraph can be supplied to the Zoltan library in one of two
//!   compressed storage formats. Although hypergraphs are often used to
//!   represent the structure of sparse matrices, the Zoltan/PHG terminology is
//!   purely in terms of vertices and hyperedges (not rows or columns). The two
//!   compressed formats are analogous to CRS and CCS for matrices.
//!
//!   In compressed hyperedge format (ZOLTAN_COMPRESSED_EDGE) a list of global
//!   hyperedge IDs is provided. Then a single list of the hypergraph pins, is
//!   provided. A pin is the connection between a vertex and a hyperedge
//!   (corresponds to a nonzero in a sparse matrix). Pins do not have separate
//!   IDs but are rather identified by the global ID of the vertex containing
//!   the pin, and implicitly also by the hyperedge ID.
//!
//!   The other format is compressed vertex (ZOLTAN_COMPRESSED_VERTEX). In this
//!   format a list of vertex global IDs is provided. Then a list of pins
//!   ordered by vertex and then by hyperedge is provided. The pin ID in this
//!   case is the global ID of the hyperedge in which the pin appears. In both
//!   formats, an array must be provided pointing to the start in the list of
//!   pins where each hyperedge or vertex begins.
//!
//!   The purpose of this query function is to tell Zoltan in which format we
//!   will supply the hypergraph, how many vertices and hyperedges there will
//!   be, and how many pins. The actual hypergraph is supplied with a query
//!   function of the type ZOLTAN_HG_CS_FN_TYPE.
//!
//!   This query function is required by all applications using the hypergraph
//!   methods of Zoltan (unless they are using the graph-based functions with
//!   hypergraph code instead).
//! \param[inout] data Pointer to user-defined data
//! \param[inout] num_lists Upon return, the number of vertices (if using
//!   compressed vertex storage) or hyperedges (if using compressed hyperedge
//!   storage) that will be supplied to Zoltan by the MPI rank.
//! \param[inout] num_nonzeros Upon return, the number of pins (connections
//!   between vertices and hyperedges) that will be supplied to Zoltan by the
//!   MPI rank.
//! \param[inout] format Upon return, the format in which we will provide the
//!   hypergraph to Zoltan. The options are ZOLTAN_COMPRESSED_EDGE and
//!   ZOLTAN_COMPRESSED_VERTEX.
//! \param[inout] ierr Error code
//! \author J. Bakosi
//******************************************************************************
{
  HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  *num_lists = hg->numMyHEdges;
  *num_nonzeros = hg->numAllNbors;

  // We will provide a compressed hyperedge (row) format. The alternative is
  // compressed vertex (column) format: ZOLTAN_COMPRESSED_VERTEX.
  *format = ZOLTAN_COMPRESSED_EDGE;
}

static void
get_hypergraph( void *data,
                int sizeGID,
                int num_edges,
                int num_nonzeros,
                int format,
                ZOLTAN_ID_PTR edgeGID,
                int *vtxPtr,
                ZOLTAN_ID_PTR vtxGID,
                int *ierr )
//******************************************************************************
//! \brief Zoltan query function that returns a hypergraph in a compressed
//!   storage (CS) format.
//! \details The size and format of the data to be returned must have been
//!   supplied to Zoltan using a ZOLTAN_HG_SIZE_CS_FN_TYPE function.
//!
//!   When a hypergraph is distributed across multiple MPI ranks, Zoltan expects
//!   that all ranks share a consistent global numbering scheme for hyperedges
//!   and vertices. Also, no two ranks should return the same pin (matrix
//!   non-zero) in this query function. (Pin ownership is unique.)
//!
//!   This query function is required by all applications using the hypergraph
//!   methods of Zoltan (unless they are using the graph-based functions with
//!   hypergraph code instead).
//!
//! \param[inout] data Pointer to user-defined data
//! \param[in] sizeGID The number of array entries used to describe a single
//!   global ID. This value is the maximum value over all ranks of the parameter
//!   NUM_GID_ENTRIES.
//! \param[in] num_edges The number of global IDs that is expected to appear on
//!   return in edgeGID. This may correspond to either vertices or
//!   (hyper-)edges.
//! \param[in] num_nonzeros The number of pins that is expected to appear on
//!   return in vtxGID.
//! \param[in] format If format is ZOLTAN_COMPRESSED_EDGE, Zoltan expects that
//!   hyperedge global IDs will be returned in edgeGID, and that vertex global
//!   IDs will be returned in vtxGID. If it is ZOLTAN_COMPRESSED_VERTEX, then
//!   vertex global IDs are expected to be returned in edgeGID and hyperedge
//!   global IDs are expected to be returned in vtxGID.
//! \param[inout] edgeGID Upon return, a list of num_edges global IDs.
//! \param[inout] vtxPtr Upon return, this array contains num_edges integers
//!   such that the number of pins specified for hyperedge j (if format is
//!   ZOLTAN_COMPRESSED_EDGE) or vertex j (if format is
//!   ZOLTAN_COMPRESSED_VERTEX) is vtxPtr[j+1]-vtxPtr[j]. If the format is
//!   ZOLTAN_COMPRESSED_EDGE, vtxPtr[j]*sizeGID is the index into the array
//!   vtxGID where edge j's pins (vertices belonging to edge j) begin; if the
//!   format is ZOLTAN_COMPRESSED_VERTEX, vtxPtr[j]*sizeGID is the index into
//!   the array vtxGID where vertex j's pins (edges to which vertex j belongs)
//!   begin. Array indices begin at zero.
//! \param[inout] vtxGID Upon return, a list of num_pins global IDs. This is the
//!   list of the pins contained in the hyperedges or vertices listed in vtxGID.
//! \param[inout] ierr Error code
//! \author J. Bakosi
//******************************************************************************
{
  HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  if ( (num_edges != hg->numMyHEdges) ||
       (num_nonzeros != hg->numAllNbors) ||
       (format != ZOLTAN_COMPRESSED_EDGE) )
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (int i=0; i<num_edges; ++i) {
    edgeGID[i] = hg->edgeGID[i];
    vtxPtr[i] = hg->nborIndex[i];
  }

  for (int i=0; i<num_nonzeros; ++i) vtxGID[i] = hg->nborGID[i];
}

static 
std::tuple< std::pair< std::vector< std::size_t >, std::vector< std::size_t > >,
            std::pair< std::vector< std::size_t >, std::vector< std::size_t > >,
            std::size_t >
createHyperGraph( const tk::UnsMesh& graph, HGRAPH_DATA& hg )
//******************************************************************************
//  Create hypergraph data structure on MPI rank zero
//! \param[in] graph Unstructured mesh graph object reference
//! \param[inout] hg Hypergraph data structure to fill
//! \return Tuple containing at 0 elements surrounding points, at 1 points
//!   surrounding points, and at 2 number of hyperedges in graph
//! \warning This function must not be called on MPI ranks other than zero.
//! \author J. Bakosi
//******************************************************************************
{
  // Get number of points from graph. The total load is taken to be proportional
  // to the number of points of the graph which is proportional to the number of
  // unique edges in the graph.
  auto npoin = graph.size();

  // Create hypergraph data structure based on mesh graph
  hg.numMyVertices = static_cast< int >( npoin );
  hg.numMyHEdges = hg.numMyVertices;
  hg.vtxGID = (ZOLTAN_ID_PTR)
    malloc(sizeof(ZOLTAN_ID_TYPE) * static_cast<std::size_t>(hg.numMyVertices));
  hg.edgeGID = (ZOLTAN_ID_PTR)
    malloc(sizeof(ZOLTAN_ID_TYPE) * static_cast<std::size_t>(hg.numMyHEdges));
  hg.nborIndex = (int*)
    malloc(sizeof(int) * static_cast<std::size_t>(hg.numMyHEdges + 1));

  // Assign global point ids
  for (int i=0; i<hg.numMyVertices; ++i)
    hg.vtxGID[ static_cast<std::size_t>(i) ] = static_cast<ZOLTAN_ID_TYPE>(i);

  // Get tetrahedron mesh graph connectivity
  const auto& inpoel = graph.tetinpoel();

  // Generate (connectivity graph) points surrounding points of graph
  auto esup = tk::genEsup( inpoel, 4 );
  auto psup = tk::genPsup( inpoel, 4, esup );
  auto& psup1 = psup.first;
  auto& psup2 = psup.second;

  // Allocate data to store the hypergraph ids. The total number of vertices or
  // neighbors in all the hyperedges of the hypergraph, nhedge = all points
  // surrounding points + number of points, since psup does not store the
  // connection to the own point, i.e., in matrix parlance, the main-diagonal.
  // In other words, here we need the number of edges in the graph, independent
  // of direction.
  auto nhedge = psup1.size() - 1 + npoin;
  hg.numAllNbors = static_cast< int >( nhedge );
  hg.nborGID = (ZOLTAN_ID_PTR)malloc(sizeof(ZOLTAN_ID_TYPE) * nhedge);

  // Fill up hypergraph edge ids and their indices
  hg.nborIndex[0] = 0;
  for (std::size_t p=0; p<npoin; ++p) {
    hg.edgeGID[p] = static_cast< ZOLTAN_ID_TYPE >( p+1 );
    // put in own point id, i.e., main diagonal
    hg.nborGID[ hg.nborIndex[p] ] = static_cast< ZOLTAN_ID_TYPE >( p );
    int j = 1;
    // put in neighbor point ids, i.e., off-diagonals
    for (auto i=psup2[p]+1; i<=psup2[p+1]; ++i, ++j) {
      hg.nborGID[ hg.nborIndex[p] + j ] =
        static_cast< ZOLTAN_ID_TYPE >( psup1[i] );
    }
    hg.nborIndex[p+1] = hg.nborIndex[p] + j;
  }

  return std::make_tuple( esup, psup, nhedge );
}

static std::size_t
emptyHyperGraph( HGRAPH_DATA& hg )
//******************************************************************************
//  Create empty hypergraph data structures for non-zero MPI ranks
//! \param[inout] hg Hypergraph data structure to fill
//! \return Number of hyperedges in graph
//! \author J. Bakosi
//******************************************************************************
{
  hg.numMyVertices = 0;
  hg.numMyHEdges = hg.numMyVertices;
  hg.numAllNbors = 0;
  hg.nborIndex = (int*)
    malloc(sizeof(int) * static_cast<std::size_t>(hg.numMyHEdges + 1));
  hg.nborIndex[0] = 0;

  return 0;
}

static void
destroyHyperGraph( HGRAPH_DATA& hg, std::size_t nhedge )
//******************************************************************************
//  Destroy hypergraph data structure
//! \param[inout] hg Hypergraph data structure to destroy
//! \param[in] nhedge Number of hyperedges in graph
//! \author J. Bakosi
//******************************************************************************
{
  if (hg.numMyVertices > 0) free( hg.vtxGID );
  if (hg.numMyHEdges > 0) free( hg.edgeGID );
  if (hg.numMyHEdges >= 0) free( hg.nborIndex );
  if (nhedge > 0) free( hg.nborGID );
}

tuple::tagged_tuple< tag::esup,  std::pair< std::vector< std::size_t >,
                                            std::vector< std::size_t > >,
                     tag::psup,  std::pair< std::vector< std::size_t >,
                                            std::vector< std::size_t > >,
                     tag::chare, std::vector< std::size_t >,
                     tag::gid, std::vector< std::size_t > >
partitionMesh( const tk::UnsMesh& graph,
               uint64_t npart,
               const tk::Print& print )
//******************************************************************************
//  Partition mesh graph using Zoltan's hypergraph algorithm in serial
//! \param[in] graph Unstructured mesh graph object reference
//! \param[in] npart Number of desired graph partitions
//! \param[in] print Pretty printer
//! \return Tagged tuple containing elements surrounding points (at tag::esup),
//!   see tk::genEsup(), points surrounding points (at tag::psup), see
//!   tk::genPsup(), array of chare ownership IDs mapping graph points to
//!   concurrent async chares (at tag::chare), and their associated global ids
//!   (at tag::gid) so that global ids on a chare are contiguous.
//! \details This function uses Zoltan to partition the mesh graph in serial. It
//!   assumes the mesh graph only exists on MPI rank 0.
//! \author J. Bakosi
//******************************************************************************
{
  // Initialize the Zoltan library
  float ver = 0.0;
  ErrChkMPI( Zoltan_Initialize( 0, nullptr, &ver ) == ZOLTAN_OK,
             "Zoltan could not be initialized" );

  // Create Zoltan data structure
  struct Zoltan_Struct *zz;
  zz = Zoltan_Create( MPI_COMM_WORLD );
  AssertMPI( zz != nullptr, "Failed to create Zoltan data structure" );

  int peid;
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );

  // Set Zoltan parameters
  Zoltan_Set_Param( zz, "DEBUG_LEVEL", "0" );
  Zoltan_Set_Param( zz, "LB_METHOD", "HYPERGRAPH" );
  Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );
  Zoltan_Set_Param( zz, "HYPERGRAPH_PACKAGE", "PHG" );
  Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "EDGE_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "RETURN_LISTS", "PART" );
  Zoltan_Set_Param( zz, "NUM_GLOBAL_PARTS", std::to_string(npart).c_str() );

  HGRAPH_DATA hg;
  std::size_t nhedge = 0;
  std::pair< std::vector< std::size_t >, std::vector< std::size_t > > esup;
  std::pair< std::vector< std::size_t >, std::vector< std::size_t > > psup;
  if (peid == 0)  
    std::tie( esup, psup, nhedge ) = createHyperGraph( graph, hg );
  else
    nhedge = emptyHyperGraph( hg );

  // Set Zoltan query functions
  Zoltan_Set_Num_Obj_Fn( zz, get_number_of_vertices, &hg );
  Zoltan_Set_Obj_List_Fn( zz, get_vertex_list, &hg );
  Zoltan_Set_HG_Size_CS_Fn( zz, get_hypergraph_size, &hg );
  Zoltan_Set_HG_CS_Fn( zz, get_hypergraph, &hg );

  // Perform partitioning using Zoltan
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids,
                exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int rc =
    Zoltan_LB_Partition(
      zz,                // Input hypergraph data structure (remaining: output)
      &changes,          // 1 if partitioning was changed, 0 otherwise
      &numGidEntries,    // Number of integers used for a global ID
      &numLidEntries,    // Number of integers used for a local ID
      &numImport,        // Number of vertices to be sent to me
      &importGlobalGids, // Global IDs of vertices to be sent to me
      &importLocalGids,  // Local IDs of vertices to be sent to me
      &importProcs,      // Process rank for source of each incoming vertex
      &importToPart,     // New partition for each incoming vertex
      &numExport,        // Number of vertices I must send to other processes
      &exportGlobalGids, // Global IDs of the vertices I must send
      &exportLocalGids,  // Local IDs of the vertices I must send
      &exportProcs,      // Process to which I send each of the vertices
      &exportToPart );   // Partition to which each vertex will belong

  // Destructor lambda
  auto destruct = [&]() {
    destroyHyperGraph( hg, nhedge );    // destroy hypergraph data structure
    Zoltan_Destroy( &zz );              // destroy Zoltan data structure
  };

  if (rc != ZOLTAN_OK) {
    destruct();
    Throw( "Zoltan_LB_Partition failed" );
  }

  // Will return, only on MPI rank 0, array of chare IDs corresponding to the
  // ownership of all points in the mesh graph, i.e., the coloring
  std::vector< std::size_t > chare;
  for( int p=0; p<numExport; ++p )
    chare.push_back( static_cast< std::size_t >( exportToPart[p] ) );

  std::size_t nchare = 1;
  if (peid == 0) {
    auto minmax = std::minmax_element( begin(chare), end(chare) );
    nchare = *minmax.second - *minmax.first + 1;

    if (npart > nchare)
      print << "\n>>> WARNING: The number of parts returned from the graph "
               "partitioner (" + std::to_string(nchare) + ") is smaller than "
               "the number of work units computed (" + std::to_string(npart) +
               ") based on the degree of virtualization desired. This may not "
               "be a problem of itself, however, it may be an indication of a "
               "too large overdecomposition. Solution 1: decrease the "
               "virtualization to a lower value using the command-line "
               "argument '-u'. Solution 2: decrease the number processing "
               "elements (PEs) using the charmrun command-line argument '+pN' "
               "where N is the number of PEs, which implicitly increases the "
               "size (and thus decreases the number) of work units.";
    else if (npart < nchare)
      Throw( "The number of parts returned from Zoltan ("
             + std::to_string(nchare) + ") is larger than the desired number "
             "of parts (" + std::to_string(npart) + ")?" );
  }

  // Construct global ids contiguous per chare
  std::vector< std::size_t > gid;
  for (std::size_t i=0; i<nchare; ++i)
    for (std::size_t n=0; n<chare.size(); ++n)
      if (chare[n] == i) gid.push_back( n );

  if (peid == 0) {
    for (auto i : chare) std::cout << i << " ";
    std::cout << '\n';
    for (auto i : gid) std::cout << i << " ";
    std::cout << '\n';
  }

  // Free the arrays allocated by Zoltan_LB_Partition
  Zoltan_LB_Free_Part( &importGlobalGids, &importLocalGids, &importProcs,
                       &importToPart );
  Zoltan_LB_Free_Part( &exportGlobalGids, &exportLocalGids, &exportProcs,
                       &exportToPart);

  // Free hypergraph and Zoltan data structure
  destruct();

  return { esup, psup, chare, gid };
}

} // zoltan::
} // tk::
