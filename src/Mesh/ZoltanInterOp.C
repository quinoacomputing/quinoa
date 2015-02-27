//******************************************************************************
/*!
  \file      src/Mesh/ZoltanInterOp.C
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:20:29 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Zoltan library
  \details   Interoperation with the Zoltan library, used for static mesh
    partitioning.
*/
//******************************************************************************

#include <zoltan.h>

#include <ZoltanInterOp.h>

namespace tk {
namespace zoltan {

void partitionMesh( const tk::UnsMesh& mesh ) {

  // Initialize the Zoltan library
  float ver = 0.0;
  ErrChk( Zoltan_Initialize( 0, nullptr, &ver ) == ZOLTAN_OK,
          "Zoltan could not be initialized" );

  // Create Zoltan data structure
  struct Zoltan_Struct *z;
  z = Zoltan_Create( MPI_COMM_WORLD );
  Assert( z != nullptr, "Zoltan_Create failed" );

  std::cout << "Will now partition...\n";

  // Destroy Zoltan data structure
  Zoltan_Destroy( &z );
}

} // zoltan::
} // tk::
