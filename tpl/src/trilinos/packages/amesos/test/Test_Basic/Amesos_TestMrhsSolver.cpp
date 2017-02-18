// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Amesos_ConfigDefs.h"
#include "Teuchos_ParameterList.hpp"
//#include "Trilinos_Util_ReadTriples2Epetra.h"
//#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Trilinos_Util.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_SLUS
#include "Epetra_SLU.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#ifdef HAVE_AMESOS_PARDISO
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
#include "Amesos_Paraklete.h"
#endif
#if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI)
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_SLUD
#include "SuperludistOO.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SLUD2
#include "Superludist2_OO.h"
#endif
#ifdef TEST_SPOOLES
#include "SpoolesOO.h"
#endif
#include "SparseSolverResult.h"
#include "Amesos_TestSolver.h"
#include "CrsMatrixTranspose.h"
#include "SparseDirectTimingVars.h"

#include <vector>
//
//  TestMrhsSolver.cpp reads in a matrix in Harwell-Boeing format, 
//  calls one of the sparse direct solvers, using multiple right hand sides
//  (one per solve) and computes the error and residual.  
//
//  TestSolver ignores the Harwell-Boeing right hand sides, creating
//  random right hand sides instead.  
//
//  TestMrhsSolver can test either A x = b or A^T x = b.
//  This can be a bit confusing because sparse direct solvers 
//  use compressed column storage - the transpose of Trilinos'
//  sparse row storage.
//
//  Matrices:
//    readA - Serial.  As read from the file.
//    transposeA - Serial.  The transpose of readA.
//    serialA - if (transpose) then transposeA else readA 
//    distributedA - readA distributed to all processes
//    passA - if ( distributed ) then distributedA else serialA
//
//
int Amesos_TestMrhsSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves, 
		     SparseSolverType SparseSolver, bool transpose, 
		     int special, AMESOS_MatrixType matrix_type ) {


  Comm.Barrier();

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;

  std::string FileName = matrix_file ;
  int FN_Size = FileName.size() ; 
  std::string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  std::string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );
  bool NonContiguousMap = false; 

  if ( LastFiveBytes == ".triU" ) { 
    // Call routine to read in unsymmetric Triplet matrix
    NonContiguousMap = true; 
    EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, false, Comm, readMap, readA, readx, 
						      readb, readxexact, NonContiguousMap ) );
  } else {
    if ( LastFiveBytes == ".triS" ) { 
      NonContiguousMap = true; 
      // Call routine to read in symmetric Triplet matrix
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, true, Comm, readMap, readA, readx, 
							readb, readxexact) );
    } else {
      if (  LastFourBytes == ".mtx" ) { 
	EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra( matrix_file, Comm, readMap, 
							       readA, readx, readb, readxexact) );
      } else {
	// Call routine to read in HB problem
	Trilinos_Util_ReadHb2Epetra( matrix_file, Comm, readMap, readA, readx, 
						     readb, readxexact) ;
      }
    }
  }


  Epetra_CrsMatrix transposeA(Copy, *readMap, 0);
  Epetra_CrsMatrix *serialA ; 

  if ( transpose ) {
    assert( CrsMatrixTranspose( readA, &transposeA ) == 0 ); 
    serialA = &transposeA ; 
  } else {
    serialA = readA ; 
  }

  
  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);
  Epetra_Map* map_;

  if( NonContiguousMap ) {
    //
    //  map gives us NumMyElements and MyFirstElement;
    //
    int NumGlobalElements =  readMap->NumGlobalElements();
    int NumMyElements = map.NumMyElements();
    int MyFirstElement = map.MinMyGID();
    std::vector<int> MapMap_( NumGlobalElements );
    readMap->MyGlobalElements( &MapMap_[0] ) ;
    Comm.Broadcast( &MapMap_[0], NumGlobalElements, 0 ) ; 
    map_ = new Epetra_Map( NumGlobalElements, NumMyElements, &MapMap_[MyFirstElement], 0, Comm);
  } else {
    map_ = new Epetra_Map( map ) ; 
  }


  // Create Exporter to distribute read-in matrix and vectors
  Epetra_Export exporter(*readMap, *map_);
  Epetra_CrsMatrix A(Copy, *map_, 0);

  Epetra_RowMatrix * passA = 0; 
  Epetra_MultiVector * passx = 0; 
  Epetra_MultiVector * passb = 0;
  Epetra_MultiVector * passxexact = 0;
  Epetra_MultiVector * passresid = 0;
  Epetra_MultiVector * passtmp = 0;

  Epetra_MultiVector x(*map_,numsolves);
  Epetra_MultiVector b(*map_,numsolves);
  Epetra_MultiVector xexact(*map_,numsolves);
  Epetra_MultiVector resid(*map_,numsolves);
  Epetra_MultiVector tmp(*map_,numsolves);


  Epetra_MultiVector serialx(*readMap,numsolves);
  Epetra_MultiVector serialb(*readMap,numsolves);
  Epetra_MultiVector serialxexact(*readMap,numsolves);
  Epetra_MultiVector serialresid(*readMap,numsolves);
  Epetra_MultiVector serialtmp(*readMap,numsolves);

  bool distribute_matrix = ( matrix_type == AMESOS_Distributed ) ; 
  if ( distribute_matrix ) { 
    //
    //  Initialize x, b and xexact to the values read in from the file
    //

    A.Export(*serialA, exporter, Add);
    Comm.Barrier();

    assert(A.FillComplete()==0);    
    Comm.Barrier();

    passA = &A; 
    passx = &x; 
    passb = &b;
    passxexact = &xexact;
    passresid = &resid;
    passtmp = &tmp;
  } else { 
    passA = serialA; 
    passx = &serialx; 
    passb = &serialb;
    passxexact = &serialxexact;
    passresid = &serialresid;
    passtmp = &serialtmp;
  }

  passxexact->SetSeed(131) ; 
  passxexact->Random();
  passx->SetSeed(11231) ; 
  passx->Random();

  passb->PutScalar( 0.0 );
  passA->Multiply( transpose, *passxexact, *passb ) ; 

  Epetra_MultiVector CopyB( *passb ) ;

  double Anorm = passA->NormInf() ; 
  SparseDirectTimingVars::SS_Result.Set_Anorm(Anorm) ;

  Epetra_LinearProblem Problem(  (Epetra_RowMatrix *) passA, 
				 (Epetra_MultiVector *) passx, 
				 (Epetra_MultiVector *) passb );

  double max_resid = 0.0;
  for ( int j = 0 ; j < special+1 ; j++ ) { 
    
    Epetra_Time TotalTime( Comm ) ; 
    if ( false ) { 
#ifdef TEST_UMFPACK

      unused code

    } else if ( SparseSolver == UMFPACK ) { 
      UmfpackOO umfpack( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
      
      umfpack.SetTrans( transpose ) ; 
      umfpack.Solve() ; 
#endif
#ifdef TEST_SUPERLU
    } else if ( SparseSolver == SuperLU ) { 
      SuperluserialOO superluserial ; 
      superluserial.SetUserMatrix( (Epetra_RowMatrix *) passA) ; 

      superluserial.SetPermc( SuperLU_permc ) ; 
      superluserial.SetTrans( transpose ) ; 
      superluserial.SetUseDGSSV( special == 0 ) ; 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	superluserial.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	superluserial.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	//      superluserial.SetRHS( (Epetra_MultiVector *) passb_i ; 
	superluserial.Solve() ; 
	if ( i == 0 ) {
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	} else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_SLUD
    } else if ( SparseSolver == SuperLUdist ) { 
      SuperludistOO superludist( Problem ) ; 
      superludist.SetTrans( transpose ) ; 

      bool factor = true; 
      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( superludist.Solve( factor ) ); 
	factor = false; 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_SLUD2
    } else if ( SparseSolver == SuperLUdist2 ) { 
      Superludist2_OO superludist2( Problem ) ; 
      superludist2.SetTrans( transpose ) ; 

      bool factor = true; 
      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( superludist2.Solve( factor ) ); 
	factor = false; 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_DSCPACK
    } else if ( SparseSolver == DSCPACK ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Dscpack dscpack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( dscpack.SetParameters( ParamList ) ); 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( dscpack.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_UMFPACK
    } else if ( SparseSolver == UMFPACK ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Umfpack umfpack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( umfpack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( umfpack.SetUseTranspose( transpose ) ); 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( umfpack.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_SUPERLU
    } else if ( SparseSolver == SUPERLU ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Superlu superlu( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( superlu.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( superlu.SetUseTranspose( transpose ) ); 

      EPETRA_CHK_ERR( superlu.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( superlu.NumericFactorization(  ) ); 
      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( superlu.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_SLUS
    } else if ( SparseSolver == SuperLU ) { 
      Epetra_SLU superluserial( &Problem ) ;
      
      bool factor = true; 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( superluserial.Solve( true, false, factor, 2, -1, true, transpose ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_KLU
    } else if ( SparseSolver == KLU ) { 
      Teuchos::ParameterList ParamList ;
      //      ParamList.set("OutputLevel",2);
      Amesos_Klu klu( Problem ) ; 
      // ParamList.set ("ScaleMethod", 0) ;
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( klu.SetParameters( ParamList ) ); 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( klu.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( klu.SetUseTranspose( transpose ) ); 

      EPETRA_CHK_ERR( klu.SymbolicFactorization(  ) ); 
      for ( int trials = 0 ; trials <= 1 ; trials++) {
	  EPETRA_CHK_ERR( klu.NumericFactorization(  ) ); 
	  for ( int i= 0 ; i < numsolves ; i++ ) {
	    //    set up to sovle A X[:,i] = B[:,i]
	    Epetra_Vector *passb_i = (*passb)(i) ;
	    Epetra_Vector *passx_i = (*passx)(i) ;
	    Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	    Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );

	    EPETRA_CHK_ERR( klu.Solve( ) ); 
	    if ( i == 0 ) {
	      SparseDirectTimingVars::SS_Result.Set_First_Time(
		      TotalTime.ElapsedTime() ); 
	    } else {
	      if ( i < numsolves-1 ) 
		SparseDirectTimingVars::SS_Result.Set_Middle_Time(
			TotalTime.ElapsedTime() ); 
	      else
		SparseDirectTimingVars::SS_Result.Set_Last_Time(
			TotalTime.ElapsedTime() ); 
	    }
	  }
      }
#endif
#ifdef HAVE_AMESOS_LAPACK
    } else if ( SparseSolver == LAPACK ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Lapack lapack( Problem ) ; 
      EPETRA_CHK_ERR( lapack.SetUseTranspose( transpose ) ); 

      EPETRA_CHK_ERR( lapack.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( lapack.NumericFactorization(  ) ); 
      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( lapack.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_TAUCS
    } else if ( SparseSolver == TAUCS ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Taucs taucs( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( taucs.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( taucs.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( taucs.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( taucs.NumericFactorization( ) ); 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( taucs.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_PARDISO
    } else if ( SparseSolver == PARDISO ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Pardiso pardiso( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( pardiso.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( pardiso.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( pardiso.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( pardiso.NumericFactorization( ) ); 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( pardiso.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_PARAKLETE
    } else if ( SparseSolver == PARAKLETE ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Paraklete paraklete( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( paraklete.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( paraklete.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( paraklete.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( paraklete.NumericFactorization( ) ); 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( paraklete.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI)
    } else if ( SparseSolver == MUMPS ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Mumps mumps( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( mumps.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( mumps.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( mumps.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( mumps.NumericFactorization( ) ); 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( mumps.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_SCALAPACK
    } else if ( SparseSolver == SCALAPACK ) { 
      Teuchos::ParameterList ParamList ;
      Amesos_Scalapack scalapack( Problem ) ; 
      ParamList.set( "MaxProcs", -3 );
      EPETRA_CHK_ERR( scalapack.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( scalapack.SetUseTranspose( transpose ) ); 

      EPETRA_CHK_ERR( scalapack.SymbolicFactorization( ) ); 
      EPETRA_CHK_ERR( scalapack.NumericFactorization( ) ); 
      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( scalapack.Solve( ) ); 
	if ( i == 0 ) 
	  SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 
	else { 
	  if ( i < numsolves-1 ) 
	    SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	  else
	    SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
	}

      }
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
    } else if ( SparseSolver == SUPERLUDIST ) { 
      Teuchos::ParameterList ParamList ;
      ParamList.set( "MaxProcs", -3 );
      Amesos_Superludist superludist( Problem ) ; 
      EPETRA_CHK_ERR( superludist.SetParameters( ParamList ) ); 
      EPETRA_CHK_ERR( superludist.SetUseTranspose( transpose ) ); 
      EPETRA_CHK_ERR( superludist.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( superludist.NumericFactorization(  ) ); 
      SparseDirectTimingVars::SS_Result.Set_First_Time( TotalTime.ElapsedTime() ); 

      for ( int i= 0 ; i < numsolves ; i++ ) { 
	//    set up to sovle A X[:,i] = B[:,i]
	Epetra_Vector *passb_i = (*passb)(i) ;
	Epetra_Vector *passx_i = (*passx)(i) ;
	Problem.SetLHS( dynamic_cast<Epetra_MultiVector *>(passx_i) ) ;
	Problem.SetRHS( dynamic_cast<Epetra_MultiVector *>(passb_i) );
	EPETRA_CHK_ERR( superludist.Solve( ) ); 
	if ( i < numsolves-1 ) 
	  SparseDirectTimingVars::SS_Result.Set_Middle_Time( TotalTime.ElapsedTime() ); 
	else
	  SparseDirectTimingVars::SS_Result.Set_Last_Time( TotalTime.ElapsedTime() ); 
      }
#endif
#ifdef TEST_SPOOLES
    } else if ( SparseSolver == SPOOLES ) { 
      SpoolesOO spooles( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
    
      spooles.SetTrans( transpose ) ; 
      spooles.Solve() ;
#endif 
#ifdef TEST_SPOOLESSERIAL
    } else if ( SparseSolver == SPOOLESSERIAL ) { 
      SpoolesserialOO spoolesserial( (Epetra_RowMatrix *) passA, 
				     (Epetra_MultiVector *) passx, 
				     (Epetra_MultiVector *) passb ) ; 
    
      spoolesserial.Solve() ;
#endif 
    } else { 
      SparseDirectTimingVars::log_file << "Solver not implemented yet" << std::endl ;
      std::cerr << "\n\n####################  Requested solver not available (Or not tested with multiple RHS) on this platform #####################\n" << std::endl ;
    }

    SparseDirectTimingVars::SS_Result.Set_Total_Time( TotalTime.ElapsedTime() ); 

    //
    //  Compute the error = norm(xcomp - xexact )
    //
    std::vector <double> error(numsolves) ; 
    double max_error = 0.0;
  
    passresid->Update(1.0, *passx, -1.0, *passxexact, 0.0);

    passresid->Norm2(&error[0]);
    for ( int i = 0 ; i< numsolves; i++ ) 
      if ( error[i] > max_error ) max_error = error[i] ; 
    SparseDirectTimingVars::SS_Result.Set_Error(max_error) ;

    //  passxexact->Norm2(&error[0] ) ; 
    //  passx->Norm2(&error ) ; 

    //
    //  Compute the residual = norm(Ax - b)
    //
    std::vector <double> residual(numsolves) ; 
  
    passtmp->PutScalar(0.0);
    passA->Multiply( transpose, *passx, *passtmp);
    passresid->Update(1.0, *passtmp, -1.0, *passb, 0.0); 
    //    passresid->Update(1.0, *passtmp, -1.0, CopyB, 0.0); 
    passresid->Norm2(&residual[0]);

    for ( int i = 0 ; i< numsolves; i++ ) 
      if ( residual[i] > max_resid ) max_resid = residual[i] ; 


    SparseDirectTimingVars::SS_Result.Set_Residual(max_resid) ;
    
    std::vector <double> bnorm(numsolves); 
    passb->Norm2( &bnorm[0] ) ; 
    SparseDirectTimingVars::SS_Result.Set_Bnorm(bnorm[0]) ;

    std::vector <double> xnorm(numsolves); 
    passx->Norm2( &xnorm[0] ) ; 
    SparseDirectTimingVars::SS_Result.Set_Xnorm(xnorm[0]) ;

  }
  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  delete map_;

  Comm.Barrier();
   return 0;
}
