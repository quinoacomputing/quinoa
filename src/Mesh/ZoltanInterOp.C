//******************************************************************************
/*!
  \file      src/Mesh/ZoltanInterOp.C
  \author    J. Bakosi
  \date      Tue 17 Mar 2015 07:27:50 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
//******************************************************************************
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <zoltan.h>
#include <charm++.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <ExceptionMPI.h>
#include <ZoltanInterOp.h>

namespace tk {
namespace zoltan {

//! Zoltan hypergraph data structure
struct HGRAPH_DATA {
  int numMyVertices;            //!< number of vertices that I own initially
  ZOLTAN_ID_TYPE *vtxGID;       //!< global ID of these vertices
  int numMyHEdges;              //!< number of my hyperedges
  int numAllNbors;              //!< number of vertices in my hyperedges
  ZOLTAN_ID_TYPE *edgeGID;      //!< global ID of each of my hyperedges
  int *nborIndex;               //!< index into nborGID array of edge's vertices
  //!< vertices of edge edgeGID[i] begin at nborGID[nborIndex[i]]
  ZOLTAN_ID_TYPE *nborGID;
};

static int get_number_of_vertices( void *data, int *ierr ) {
  HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;
  return hg->numMyVertices;
}

void partitionMesh( const tk::UnsMesh& mesh ) {

  // Initialize the Zoltan library
  float ver = 0.0;
  ErrChkMPI( Zoltan_Initialize( 0, nullptr, &ver ) == ZOLTAN_OK,
             "Zoltan could not be initialized" );

  // Create Zoltan data structure
  struct Zoltan_Struct *zz;
  zz = Zoltan_Create( MPI_COMM_WORLD );
  AssertMPI( zz != nullptr, "Failed to create Zoltan data structure" );

  // Set Zoltan parameters
  char global_parts[10];
  sprintf( global_parts, "%d", CkNumPes() );
  Zoltan_Set_Param( zz, "DEBUG_LEVEL", "0" );
  Zoltan_Set_Param( zz, "LB_METHOD", "HYPERGRAPH" );
  Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );
  Zoltan_Set_Param( zz, "HYPERGRAPH_PACKAGE", "PHG" );
  Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
  Zoltan_Set_Param( zz, "RETURN_LISTS", "PART" );
  Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "EDGE_WEIGHT_DIM", "0" );
  Zoltan_Set_Param( zz, "NUM_GLOBAL_PARTS", global_parts );

  // Create hypergraph data structure based on mesh
  HGRAPH_DATA hg;
  hg.numMyVertices = static_cast< int >( mesh.nnode() );
  hg.numMyHEdges = hg.numMyVertices;
  hg.vtxGID = (ZOLTAN_ID_PTR)malloc(sizeof(ZOLTAN_ID_TYPE) *
                                    static_cast<std::size_t>(hg.numMyVertices));
  hg.edgeGID = (ZOLTAN_ID_PTR)malloc(sizeof(ZOLTAN_ID_TYPE) *
                                     static_cast<std::size_t>(hg.numMyHEdges));
  hg.nborIndex = (int*)malloc(sizeof(int) *
                              static_cast<std::size_t>(hg.numMyHEdges + 1));
  for (int i=0; i<hg.numMyVertices; ++i)
    hg.vtxGID[ static_cast<std::size_t>(i) ] = static_cast<ZOLTAN_ID_TYPE>(i+1);

  // Set Zoltan query functions
  Zoltan_Set_Num_Obj_Fn( zz, get_number_of_vertices, &hg );

  // Destroy hypergraph data structure
  if (hg.numMyVertices > 0) {
    free( hg.vtxGID );
    free( hg.edgeGID );
    free( hg.nborIndex );
  }

  // Destroy Zoltan data structure
  Zoltan_Destroy( &zz );
}

} // zoltan::
} // tk::
