/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include <stdio.h>
#include <assert.h>

#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"

#include "CEpetra_SerialComm.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"

void stall( )
{
  long int j, stall = 32*1024*1024;
  double dummy1 = 1.0;
  double dummy2 = 1.01;

  for (j=0; j<stall; j++) {
    /* Stall so that the human can pull the plug before we
     * go up in smoke, if we are leaking */
    dummy1 *= dummy2;
    dummy2 *= dummy1;
  }

  return;
}

int main( int argc, char* argv[] )
{
  int numGlobalElements, numMyElements, indexBase;
  int numGlobalElements_rtn, numIterations, i;

  numGlobalElements = 64*128*1024;

  if (argc > 1) {
    numIterations = atoi(argv[1]);
  } else {
    printf( "Syntax: %s <numIters>\n", argv[0] );
    printf( "\nEnd Result: TEST FAILED\n" );
    return 1;
  }

  CT_Epetra_SerialComm_ID_t scommID;
  CT_Epetra_Comm_ID_t commID;

  CT_Epetra_Map_ID_t mapID;
  CT_Epetra_BlockMap_ID_t bmapID;

  CT_Epetra_Vector_ID_t xID, bID;

  printf( "numIterations = %d\n", numIterations);

  /* Create an Epetra_SerialComm and cast to an Epetra_Comm so that
   * it can be passed to functions expecting the latter */
  scommID = Epetra_SerialComm_Create();
  commID = Epetra_Comm_Cast(Epetra_SerialComm_Abstract(scommID));

  /* Create an Epetra_Map and cast to an Epetra_BlockMap so that
   * a) it can be passed to functions expecting the latter and
   * b) methods implemented only in BlockMap can be invoked on the Map */
  indexBase = 0;  /* use indexBase = 0 unless you know what you're doing! */
  mapID = Epetra_Map_Create(numGlobalElements, indexBase, commID);
  bmapID = Epetra_BlockMap_Cast(Epetra_Map_Abstract(mapID));

  /* Check the properties of the map */
  numGlobalElements_rtn = Epetra_BlockMap_NumGlobalElements(bmapID);
  printf( "NumGlobalElements = %d\n", numGlobalElements_rtn );
  assert( numGlobalElements == numGlobalElements_rtn );

  numMyElements = Epetra_BlockMap_NumMyElements(bmapID);
  printf( "NumMyElements = %d\n", numMyElements );

  printf("\n");
  for (i=0; i<numIterations; i++) {
    printf( "Iteration %d...\n", i );

    /* Create an Epetra_Vector and then destroy it */
    /* If the tables leak memory, we humans should notice */
    xID = Epetra_Vector_Create(bmapID, TRUE);  /* zero it too */
    stall();  stall();  stall();
    Epetra_Vector_Destroy(&xID);
    stall();
  }

  printf( "\nEnd Result: TEST PASSED\n" );

  return 0;
}
