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


#include <iostream>

#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"
#include "CEpetra_Comm.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm.h"
#endif
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Import.h"
#include "CEpetra_Export.h"


int checkresult(int MyResult, int Correct, CT_Epetra_Comm_ID_t CommID, int MyPID, const char *pc)
{
  int Failed = (MyResult == Correct ? 0 : 1);

  if (Failed != 0)
    std::cout << "  " << pc << " FAILED on processor " << MyPID << ":  " << MyResult << " != " << Correct << std::endl;

  int SomeFailed;
  Epetra_Comm_MaxAll_Int(CommID, &Failed, &SomeFailed, 1);

  if ((SomeFailed == 0) && (MyPID == 0))
    std::cout << "  " << pc << " passed!" << std::endl;
  Epetra_Comm_Barrier(CommID);

  return SomeFailed;
}

int checkresultarray(int MyResult[], int Correct[], int Count, CT_Epetra_Comm_ID_t CommID, int MyPID, const char *pc)
{
  int Failed = 0;
  for (int i=0; i<Count; i++) {
    if (MyResult[i] != Correct[i]) Failed = 1;
  }

  if (Failed != 0) {
    std::cout << "  " << pc << " FAILED on processor " << MyPID << ": {" << MyResult[0];
    for (int i=1; i<Count; i++) std::cout << "," << MyResult[i];
    std::cout << "} != {" << Correct[0];
    for (int i=1; i<Count; i++) std::cout << "," << Correct[i];
    std::cout << "}" << std::endl;
  }

  int SomeFailed;
  Epetra_Comm_MaxAll_Int(CommID, &Failed, &SomeFailed, 1);

  if ((SomeFailed == 0) && (MyPID == 0))
    std::cout << "  " << pc << " passed!" << std::endl;
  Epetra_Comm_Barrier(CommID);

  return SomeFailed;
}

int main(int argc, char *argv[])
{
  int success = 0;

#ifdef HAVE_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);

  const int ForceNumProc = 3;

  /* Set up communication */
  CT_Epetra_Comm_ID_t CommID = Epetra_Comm_Degeneralize(Epetra_MpiComm_Generalize(
      Epetra_MpiComm_Create( MPI_COMM_WORLD )));

  int MyPID = Epetra_Comm_MyPID(CommID);
  int NumProc = Epetra_Comm_NumProc(CommID);

  if (NumProc != ForceNumProc) {
    std::cerr << "Run this with " << ForceNumProc << " processors!" << std::endl;
    std::cout << "End Result: TEST FAILED" << std::endl;
    success = 1;
    return success;
  }

  Epetra_Comm_Barrier(CommID);

  /* Create the source map */
  if (MyPID == 0) std::cout << "Creating source map..." << std::endl;
  int IndexBase = 0;
  const int NumMyElements = 3;
  int NumGlobalElements = NumMyElements * NumProc;
  int off = NumMyElements*MyPID;
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
  CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
      Epetra_Map_Create_Arbitrary(NumGlobalElements, NumMyElements,
          MyGlobalElements, IndexBase, CommID)));

  Epetra_Comm_Barrier(CommID);

  /* Create the target map */
  if (MyPID == 0) std::cout << "Creating target map..." << std::endl;
  const int NumMyElements2 = 5;
  int NumGlobalElements2 = NumMyElements2 * NumProc;
  int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
  int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
  int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
  int MyGlobalElements2[NumMyElements2];
  for (int i=0; i<NumMyElements2; i++) {
    switch (MyPID) {
    case 0:
      MyGlobalElements2[i] = MyGlobalElements2a[i];
      break;
    case 1:
      MyGlobalElements2[i] = MyGlobalElements2b[i];
      break;
    case 2:
      MyGlobalElements2[i] = MyGlobalElements2c[i];
      break;
    }
  }
  CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Degeneralize(Epetra_Map_Generalize(
       Epetra_Map_Create_Arbitrary(NumGlobalElements2, NumMyElements2,
           MyGlobalElements2, IndexBase, CommID)));

  Epetra_Comm_Barrier(CommID);

{ /* IMPORT BLOCK */

  /* Create an importer */
  if (MyPID == 0) std::cout << "Creating importer..." << std::endl;
  CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(btarID, bsrcID);

  Epetra_Comm_Barrier(CommID);

  /* Try out the wrapper functions... */
  if (MyPID == 0) std::cout << "Testing importer..." << std::endl;

  Epetra_Comm_Barrier(CommID);

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumSameIDs(selfID);
    int Correct[ForceNumProc] = {3, 0, 0};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSameIDs");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumPermuteIDs(selfID);
    int Correct[ForceNumProc] = {0, 3, 3};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumPermuteIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_PermuteToLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {1, 2, 3};
      int Correct2c[MaxCount] = {2, 3, 4};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteToLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_PermuteFromLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {0, 1, 2};
      int Correct2c[MaxCount] = {0, 1, 2};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteFromLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumRemoteIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRemoteIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_RemoteLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {3, 4};
      int Correct2b[MaxCount] = {0, 4};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "RemoteLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumExportIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumExportIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_ExportLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {2, 0};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {2, 0};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Import_ExportPIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {1, 2};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportPIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumSend(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSend");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Import_NumRecv(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRecv");
  }

  /****************************************************************/

  Epetra_Comm_Barrier(CommID);

  /* Destroy importer */
  Epetra_Import_Destroy(&selfID);

  Epetra_Comm_Barrier(CommID);

} /* IMPORT BLOCK */

{ /* EXPORT BLOCK */

  /* Create an exporter */
  if (MyPID == 0) std::cout << "Creating exporter..." << std::endl;
  /* Intentionally swapping source and target maps! */
  CT_Epetra_Export_ID_t selfID = Epetra_Export_Create(btarID, bsrcID);

  Epetra_Comm_Barrier(CommID);

  /* Try out the wrapper functions... */
  if (MyPID == 0) std::cout << "Testing exporter..." << std::endl;

  Epetra_Comm_Barrier(CommID);

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumSameIDs(selfID);
    int Correct[ForceNumProc] = {3, 0, 0};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSameIDs");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumPermuteIDs(selfID);
    int Correct[ForceNumProc] = {0, 3, 3};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumPermuteIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_PermuteToLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {0, 1, 2};
      int Correct2c[MaxCount] = {0, 1, 2};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteToLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_PermuteFromLIDs(selfID);
      const int MaxCount = 3;
      int Correct2a[MaxCount];
      int Correct2b[MaxCount] = {1, 2, 3};
      int Correct2c[MaxCount] = {2, 3, 4};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "PermuteFromLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumRemoteIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRemoteIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_RemoteLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {2, 0};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {2, 0};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "RemoteLIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumExportIDs(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    int success1 = checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumExportIDs");
    success += success1;

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_ExportLIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {3, 4};
      int Correct2b[MaxCount] = {0, 4};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportLIDs");
    }

    if (success1 == 0) {
      int *MyResult2 = Epetra_Export_ExportPIDs(selfID);
      const int MaxCount = 2;
      int Correct2a[MaxCount] = {1, 2};
      int Correct2b[MaxCount] = {0, 2};
      int Correct2c[MaxCount] = {0, 1};
      int *Correct2 = (MyPID == 0 ? Correct2a : (MyPID == 1 ? Correct2b : Correct2c));
      success += checkresultarray(MyResult2, Correct2, MyResult, CommID, MyPID, "ExportPIDs");
    }
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumSend(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumSend");
  }

  /****************************************************************/

  {
    int MyResult = Epetra_Export_NumRecv(selfID);
    int Correct[ForceNumProc] = {2, 2, 2};
    success += checkresult(MyResult, Correct[MyPID], CommID, MyPID, "NumRecv");
  }

  /****************************************************************/

  Epetra_Comm_Barrier(CommID);

  /* Destroy exporter */
  Epetra_Export_Destroy(&selfID);

  Epetra_Comm_Barrier(CommID);

} /* EXPORT BLOCK */

  if (success == 0)
    std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
  else
    std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
  
  MPI_Finalize() ;

#else /* HAVE_MPI */

  std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
  return 1;

#endif /* HAVE_MPI */

  return ((success == 0) ? 0 : 1);
}

