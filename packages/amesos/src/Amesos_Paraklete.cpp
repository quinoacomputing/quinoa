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

#include "Amesos_Paraklete.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Amesos_Support.h"
extern "C" {
#include "amesos_paraklete_decl.h"
}

// The following hack allows assert to work even though it has been turned off by paraklete.h bug #1952 

#undef	_ASSERT_H
#undef	assert
#undef __ASSERT_VOID_CAST
#undef NDEBUG
#include <assert.h>

namespace {

using Teuchos::RCP;

template<class T, class DeleteFunctor>
class DeallocFunctorDeleteWithCommon
{
public:
  DeallocFunctorDeleteWithCommon(
				 const RCP<paraklete_common>  &common
				 ,DeleteFunctor                        deleteFunctor
				 )
    : common_(common), deleteFunctor_(deleteFunctor)
    {}
  typedef T ptr_t;
  void free( T* ptr ) {
    if(ptr) deleteFunctor_(&ptr,&*common_);
  }
private:
  Teuchos::RCP<paraklete_common> common_;
  DeleteFunctor deleteFunctor_;
  DeallocFunctorDeleteWithCommon(); // Not defined and not to be called!
};

template<class T, class DeleteFunctor>
DeallocFunctorDeleteWithCommon<T,DeleteFunctor>
deallocFunctorDeleteWithCommon(
			       const RCP<paraklete_common>  &common
			       ,DeleteFunctor                        deleteFunctor
			       )
{
  return DeallocFunctorDeleteWithCommon<T,DeleteFunctor>(common,deleteFunctor);
}


} // namespace 



class Amesos_Paraklete_Pimpl {
public:


  cholmod_sparse pk_A_ ;
  Teuchos::RCP<paraklete_symbolic> LUsymbolic_ ;
  Teuchos::RCP<paraklete_numeric> LUnumeric_ ;
  Teuchos::RCP<paraklete_common> common_;



} ;

//=============================================================================
Amesos_Paraklete::Amesos_Paraklete(const Epetra_LinearProblem &prob ) :
  PrivateParakleteData_( rcp( new Amesos_Paraklete_Pimpl() ) ),
  ParakleteComm_(0),
  CrsMatrixA_(0),
  TrustMe_(false),
  UseTranspose_(false),
  Problem_(&prob),
  MtxConvTime_(-1),
  MtxRedistTime_(-1),
  VecRedistTime_(-1),
  SymFactTime_(-1),
  NumFactTime_(-1),
  SolveTime_(-1),
  OverheadTime_(-1)
{

  // MS // move declaration of Problem_ above because I need it
  // MS // set up before calling Comm()
  MaxProcesses_ = -3;    // Use all processes unless the user requests otherwise
  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ;
}

//=============================================================================
Amesos_Paraklete::~Amesos_Paraklete(void) {


    if(ParakleteComm_)  {
      MPI_Comm_free(&ParakleteComm_);
      ParakleteComm_ = 0 ; 
    }

  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();
}

//=============================================================================
int Amesos_Paraklete::ExportToSerial() 
{
  if (  AddZeroToDiag_==0 && numentries_ != RowMatrixA_->NumGlobalNonzeros()) { 
    std::cerr << " The number of non zero entries in the matrix has changed since the last call to SymbolicFactorization().  " ;
    AMESOS_CHK_ERR( -2 );
  }
  if ( UseDataInPlace_ == 0 ) {
    assert ( RowMatrixA_ != 0 ) ; 
    if (  AddZeroToDiag_==0 && numentries_ != RowMatrixA_->NumGlobalNonzeros()) { 
      std::cerr << " The number of non zero entries in the matrix has changed since the last call to SymbolicFactorization().  " ;
      AMESOS_CHK_ERR( -2 );
    }

#ifdef HAVE_AMESOS_EPETRAEXT
    if ( transposer_.get() != 0  ) {
      //      int OriginalTracebackMode = CrsMatrixA_->GetTracebackMode();
      //      CrsMatrixA_->SetTracebackMode( EPETRA_MIN( OriginalTracebackMode, 0) );
      transposer_->fwd() ;
      //      CrsMatrixA_->SetTracebackMode( OriginalTracebackMode );
    }
#endif
    assert ( ImportToSerial_.get() != 0 ) ; 
    AMESOS_CHK_ERR(SerialCrsMatrixA_->Import(*StdIndexMatrix_, 
					     *ImportToSerial_, Insert ));

    if (AddZeroToDiag_ ) { 
      int OriginalTracebackMode = SerialCrsMatrixA_->GetTracebackMode() ; 
      SerialCrsMatrixA_->SetTracebackMode( EPETRA_MIN( OriginalTracebackMode, 0) ) ; // ExportToSerial is called both by PerformSymbolicFactorization() and PerformNumericFactorization().  When called by the latter, the call to insertglobalvalues is both unnecessary and illegal.  Fortunately, Epetra allows us to ignore the error message by setting the traceback mode to 0.

      //
      //  Add 0.0 to each diagonal entry to avoid empty diagonal entries;
      //
      double zero = 0.0;
      for ( int i = 0 ; i < SerialMap_->NumGlobalElements(); i++ ) 
	if ( SerialCrsMatrixA_->LRID(i) >= 0 ) 
	  SerialCrsMatrixA_->InsertGlobalValues( i, 1, &zero, &i ) ;
      SerialCrsMatrixA_->SetTracebackMode( OriginalTracebackMode ) ; 
    }
    AMESOS_CHK_ERR(SerialCrsMatrixA_->FillComplete());
    //AMESOS_CHK_ERR(SerialCrsMatrixA_->OptimizeStorage());

    if( !AddZeroToDiag_ && numentries_ != SerialMatrix_->NumGlobalNonzeros()) {
      std::cerr << " Amesos_Paraklete cannot handle this matrix.  " ;
      if ( Reindex_ ) {
	std::cerr << "Unknown error" << std::endl ; 
	AMESOS_CHK_ERR( -5 );
      } else {
	std::cerr << " Try setting the Reindex parameter to true. " << std::endl ; 
#ifndef HAVE_AMESOS_EPETRAEXT
	std::cerr << " You will need to rebuild the Amesos library with the EpetraExt library to use the reindex feature" << std::endl ; 
	std::cerr << " To rebuild Amesos with EpetraExt, add --enable-epetraext to your configure invocation" << std::endl ;
#endif
	AMESOS_CHK_ERR( -3 );
      }
    }
  }
  return 0;
}
//=============================================================================
//
// CreateLocalMatrixAndExporters() is called only by SymbolicFactorization() 
// for CrsMatrix objects.  All objects should be created here.  No assumptions 
// are made about the input operator.  I.e. it can be completely different from 
// the matrix at the time of the previous call to SymbolicFactorization(). 
//

int Amesos_Paraklete::CreateLocalMatrixAndExporters() 
{
  ResetTimer(0);

  RowMatrixA_ = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  if (RowMatrixA_ == 0) AMESOS_CHK_ERR(-1);
  
  const Epetra_Map &OriginalMatrixMap = RowMatrixA_->RowMatrixRowMap() ;
  const Epetra_Map &OriginalDomainMap = 
    UseTranspose()?GetProblem()->GetOperator()->OperatorRangeMap():
    GetProblem()->GetOperator()->OperatorDomainMap();
  const Epetra_Map &OriginalRangeMap = 
    UseTranspose()?GetProblem()->GetOperator()->OperatorDomainMap():
    GetProblem()->GetOperator()->OperatorRangeMap();
  
  NumGlobalElements_ = RowMatrixA_->NumGlobalRows();
  numentries_ = RowMatrixA_->NumGlobalNonzeros();
  assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );
  
  //
  //  Create a serial matrix
  //
  assert( NumGlobalElements_ == OriginalMatrixMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (MyPID_==0) NumMyElements_ = NumGlobalElements_;
  
  //
  //  UseDataInPlace_ is set to 1 (true) only if everything is perfectly
  //  normal.  Anything out of the ordinary reverts to the more expensive
  //  path. 
  //
  UseDataInPlace_ = ( OriginalMatrixMap.NumMyElements() ==
		      OriginalMatrixMap.NumGlobalElements() )?1:0;
  if ( ! OriginalRangeMap.SameAs( OriginalMatrixMap ) ) UseDataInPlace_ = 0 ; 
  if ( ! OriginalDomainMap.SameAs( OriginalMatrixMap ) ) UseDataInPlace_ = 0 ; 
  if ( AddZeroToDiag_ ) UseDataInPlace_ = 0 ; 
  Comm().Broadcast( &UseDataInPlace_, 1, 0 ) ;
  
  UseDataInPlace_ = 0 ; // bug - remove this someday.
    
  //
  //  Reindex matrix if necessary (and possible - i.e. CrsMatrix)
  //
  //  For now, since I don't know how to determine if we need to reindex the matrix, 
  //  we allow the user to control re-indexing through the "Reindex" parameter.
  //
  CrsMatrixA_ = dynamic_cast<Epetra_CrsMatrix *>(Problem_->GetOperator());
  Epetra_CrsMatrix* CcsMatrixA = 0 ; 
    
  //
  //  CcsMatrixA points to a matrix in Compressed Column Format 
  //  i.e. the format needed by Paraklete.  If we are solving  
  //  A^T x = b, CcsMatrixA = CrsMatrixA_, otherwise we must 
  //  transpose the matrix.
  //
  if (UseTranspose()) {
    CcsMatrixA = CrsMatrixA_ ;
  } else {
    if ( CrsMatrixA_ == 0 ) AMESOS_CHK_ERR( -7 );   //  Amesos_Paraklete only supports CrsMatrix objects in the non-transpose case
#ifdef HAVE_AMESOS_EPETRAEXT
    bool MakeDataContiguous = true;
    transposer_ = rcp ( new EpetraExt::RowMatrix_Transpose( MakeDataContiguous ));
    
    int OriginalTracebackMode = CrsMatrixA_->GetTracebackMode();
    CrsMatrixA_->SetTracebackMode( EPETRA_MIN( OriginalTracebackMode, 0) );
    
    CcsMatrixA = &(dynamic_cast<Epetra_CrsMatrix&>(((*transposer_)(*CrsMatrixA_))));
    CrsMatrixA_->SetTracebackMode( OriginalTracebackMode );
#else
    std::cerr << "Amesos_Paraklete requires the EpetraExt library to solve non-transposed problems. " << std::endl ; 
    std::cerr << " To rebuild Amesos with EpetraExt, add --enable-epetraext to your configure invocation" << std::endl ;
    AMESOS_CHK_ERR( -3 );
#endif
  }
  
  
  if  ( Reindex_ ) {
    if ( CcsMatrixA == 0 ) AMESOS_CHK_ERR(-4);
#ifdef HAVE_AMESOS_EPETRAEXT
    const Epetra_Map& OriginalMap = CcsMatrixA->RowMap();
    StdIndex_ = rcp( new Amesos_StandardIndex( OriginalMap  ) );
    //const Epetra_Map& OriginalColMap = CcsMatrixA->RowMap();
    StdIndexDomain_ = rcp( new Amesos_StandardIndex( OriginalDomainMap  ) );
    StdIndexRange_ = rcp( new Amesos_StandardIndex( OriginalRangeMap  ) );
    
    StdIndexMatrix_ = StdIndex_->StandardizeIndex( CcsMatrixA );
#else
    std::cerr << "Amesos_Paraklete requires EpetraExt to reindex matrices." << std::endl 
	 <<  " Please rebuild with the EpetraExt library by adding --enable-epetraext to your configure invocation" << std::endl ;
    AMESOS_CHK_ERR(-4);
#endif
  } else { 
    if ( UseTranspose() ){ 
      StdIndexMatrix_ = RowMatrixA_ ;
    } else {
      StdIndexMatrix_ = CcsMatrixA ;
    } 
  }
  
  //
  //  At this point, StdIndexMatrix_ points to a matrix with 
  //  standard indexing.  
  //

  //
  //  Convert Original Matrix to Serial (if it is not already)
  //
  if (UseDataInPlace_ == 1) {
    SerialMatrix_ = StdIndexMatrix_;
  } else {
    SerialMap_ = rcp(new Epetra_Map(NumGlobalElements_, NumMyElements_, 0, Comm()));
    
    if ( Problem_->GetRHS() ) 
      NumVectors_ = Problem_->GetRHS()->NumVectors() ; 
    else
      NumVectors_ = 1 ; 
    SerialXextract_ = rcp( new Epetra_MultiVector(*SerialMap_,NumVectors_));
    SerialBextract_ = rcp (new Epetra_MultiVector(*SerialMap_,NumVectors_));

    ImportToSerial_ = rcp(new Epetra_Import ( *SerialMap_, StdIndexMatrix_->RowMatrixRowMap() ) );
    if (ImportToSerial_.get() == 0) AMESOS_CHK_ERR(-5);

    SerialCrsMatrixA_ = rcp( new Epetra_CrsMatrix(Copy, *SerialMap_, 0) ) ;
    SerialMatrix_ = &*SerialCrsMatrixA_ ;
  }
  
  MtxRedistTime_ = AddTime("Total matrix redistribution time", MtxRedistTime_, 0);
  
  return(0);
}

//=============================================================================
//
//  See also pre and post conditions in Amesos_Paraklete.h
//  Preconditions:
//    firsttime specifies that this is the first time that 
//    ConertToParakleteCrs has been called, i.e. in symbolic factorization.  
//    No data allocation should happen unless firsttime=true.
//    SerialMatrix_ points to the matrix to be factored and solved
//    NumGlobalElements_ has been set to the dimension of the matrix
//    numentries_ has been set to the number of non-zeros in the matrix
//      (i.e. CreateLocalMatrixAndExporters() has been callded)
//
//  Postconditions:
//    pk_A_ contains the matrix as Paraklete needs it
//
//
int Amesos_Paraklete::ConvertToParakleteCRS(bool firsttime)
{
  ResetTimer(0);
  
  //
  //  Convert matrix to the form that Klu expects (Ap, Ai, VecAval)
  //

  if (MyPID_==0) {
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalRows());
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalCols());

    if ( ! AddZeroToDiag_ ) {
      assert( numentries_ == SerialMatrix_->NumGlobalNonzeros()) ;
    } else {
      numentries_ = SerialMatrix_->NumGlobalNonzeros() ;
    }

    Epetra_CrsMatrix *CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(SerialMatrix_);
    bool StorageOptimized = ( CrsMatrix != 0 && CrsMatrix->StorageOptimized() );

    if ( AddToDiag_ != 0.0 ) StorageOptimized = false ;

    if ( firsttime ) { 
      Ap.resize( NumGlobalElements_+1 );
      Ai.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
      if ( ! StorageOptimized ) { 
	VecAval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
	Aval = &VecAval[0];
      }
    }

    double *RowValues;
    int *ColIndices;
    int NumEntriesThisRow;

    if( StorageOptimized ) {
      if ( firsttime ) {
	Ap[0] = 0;
	for ( int MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
	  if( CrsMatrix->
	      ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues,
				ColIndices ) != 0 ) 
	    AMESOS_CHK_ERR( -6 );
          for ( int j=0; j < NumEntriesThisRow; j++ ) {
            Ai[Ap[MyRow]+j] = ColIndices[j];
          }
	  if ( MyRow == 0 ) {
	    Aval = RowValues ;
	  }
	  Ap[MyRow+1] = Ap[MyRow] + NumEntriesThisRow ;
	}
      }
    } else { 
      
      int Ai_index = 0 ;
      int MyRow;
      
      int MaxNumEntries_ = SerialMatrix_->MaxNumEntries();
      if ( firsttime && CrsMatrix == 0 ) {
	ColIndicesV_.resize(MaxNumEntries_);
	RowValuesV_.resize(MaxNumEntries_);
      }
      
      for ( MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
	if ( CrsMatrix != 0 ) {
	  if( CrsMatrix->
	      ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues,
				ColIndices ) != 0 ) 
	    AMESOS_CHK_ERR( -6 );

	} else {
	  if( SerialMatrix_->
			  ExtractMyRowCopy( MyRow, MaxNumEntries_,
					    NumEntriesThisRow, &RowValuesV_[0],
					    &ColIndicesV_[0] ) != 0 ) 
	    AMESOS_CHK_ERR( -6 );

	  RowValues =  &RowValuesV_[0];
	  ColIndices = &ColIndicesV_[0];
	}
	
	if ( firsttime ) {
	  Ap[MyRow] = Ai_index ;
	  for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	    Ai[Ai_index] = ColIndices[j] ;
	    //  VecAval[Ai_index] = RowValues[j] ;      //  We have to do this because of the hacks to get around bug #1502 
	    if (ColIndices[j] == MyRow) {
	      VecAval[Ai_index] += AddToDiag_;    
	    }
	    Ai_index++;
	  }
	} else { 
	  for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	    VecAval[Ai_index] = RowValues[j] ;     
	    if (ColIndices[j] == MyRow) {
	      VecAval[Ai_index] += AddToDiag_;   
	    }
	    Ai_index++;
	  }
	}
      }
      Ap[MyRow] = Ai_index ;
    }
  PrivateParakleteData_->pk_A_.nrow = NumGlobalElements_ ; 
  cholmod_sparse& pk_A =  PrivateParakleteData_->pk_A_ ;
  pk_A.nrow = NumGlobalElements_ ; 
  pk_A.ncol = NumGlobalElements_ ; 
  pk_A.nzmax = numentries_ ; 
  pk_A.p = &Ap[0] ;
  pk_A.i = &Ai[0] ;
  pk_A.nz = 0; 
  if ( firsttime ) { 
    pk_A.x = 0 ; 
    pk_A.xtype = CHOLMOD_PATTERN ; 
  }
  else
  {
    pk_A.x = Aval ; 
    pk_A.xtype = CHOLMOD_REAL ; 
  }

  pk_A.z = 0 ; 
  pk_A.stype = 0 ; //  symmetric 
  pk_A.itype = CHOLMOD_LONG ; 
  pk_A.dtype = CHOLMOD_DOUBLE ; 
  pk_A.sorted = 0 ; 
  pk_A.packed = 1 ; 
  }

  MtxConvTime_ = AddTime("Total matrix conversion time", MtxConvTime_, 0);

  return 0;
}

//=============================================================================
int Amesos_Paraklete::SetParameters( Teuchos::ParameterList &ParameterList ) {
  
  // ========================================= //
  // retrive PARAKLETE's parameters from list.       //
  // default values defined in the constructor //
  // ========================================= //
  
  // retrive general parameters
  
  SetStatusParameters( ParameterList );
  SetControlParameters( ParameterList );
  
  
  
  if (ParameterList.isParameter("TrustMe") ) 
    TrustMe_ = ParameterList.get<bool>( "TrustMe" );
  
#if 0
  
  unused for now
    
    if (ParameterList.isSublist("Paraklete") ) {
      Teuchos::ParameterList ParakleteParams = ParameterList.sublist("Paraklete") ;
    }
#endif
  
  return 0;
}


//=============================================================================
int Amesos_Paraklete::PerformSymbolicFactorization() 
{
  ResetTimer(0);
  
  if ( IamInGroup_ ) 
    PrivateParakleteData_->LUsymbolic_ =
      rcp( amesos_paraklete_analyze ( &PrivateParakleteData_->pk_A_, &*PrivateParakleteData_->common_ )
	   ,deallocFunctorDeleteWithCommon<paraklete_symbolic>(PrivateParakleteData_->common_,
							       amesos_paraklete_free_symbolic)
	   ,true
	   );
  
  SymFactTime_ = AddTime("Total symbolic factorization time", SymFactTime_, 0);
  
  return 0;
}

//=============================================================================
int Amesos_Paraklete::PerformNumericFactorization( ) 
{
  // Changed this; it was "if (!TrustMe)...
  // The behavior is not intuitive. Maybe we should introduce a new
  // parameter, FastSolvers or something like that, that does not perform
  // any AddTime, ResetTimer, GetTime.
  
  ResetTimer(0);
  
  //bool factor_with_pivoting = true ;
  
  if ( IamInGroup_ ) 
    PrivateParakleteData_->LUnumeric_ =
      rcp( amesos_paraklete_factorize ( &PrivateParakleteData_->pk_A_,
                                 &*PrivateParakleteData_->LUsymbolic_, 
                                 &*PrivateParakleteData_->common_ ) 
           ,deallocFunctorDeleteWithCommon<paraklete_numeric>(PrivateParakleteData_->common_,amesos_paraklete_free_numeric)
	   ,true
	   );

  NumFactTime_ = AddTime("Total numeric factorization time", NumFactTime_, 0);
  
  return 0;
}

//=============================================================================
bool Amesos_Paraklete::MatrixShapeOK() const {
  bool OK = true;
  
  // Comment by Tim:  The following code seems suspect.  The variable "OK"
  // is not set if the condition is true.
  // Does the variable "OK" default to true?
  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() !=
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) {
    OK = false;
  }
  return OK;
}


extern "C" {
  void my_handler (int status, char *file, int line, char *msg)
  {
    printf ("Error handler: %s %d %d: %s\n", file, line, status, msg) ;
    if (status != CHOLMOD_OK)
      {
	fprintf (stderr, "\n\n*********************************************\n");
	fprintf (stderr, "**** Test failure: %s %d %d %s\n", file, line,
		 status, msg) ;
	fprintf (stderr, "*********************************************\n\n");
	fflush (stderr) ;
	fflush (stdout) ;
      }
  }
}  

//=============================================================================
int Amesos_Paraklete::SymbolicFactorization() 
{
  MyPID_    = Comm().MyPID();
  NumProcs_ = Comm().NumProc();
  
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;
  
#ifdef HAVE_AMESOS_EPETRAEXT
  transposer_ = static_cast<Teuchos::ENull>( 0 ); 
#endif
  
  CreateTimer(Comm(), 2);
  
  ResetTimer(1);
  
  // "overhead" time for the following method is considered here
  AMESOS_CHK_ERR( CreateLocalMatrixAndExporters() ) ;
  assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );
  
  
  SetMaxProcesses(MaxProcesses_, *RowMatrixA_);
  
  //
  //  Perform checks in SymbolicFactorization(), but none in 
  //  NumericFactorization() or Solve()
  //
  assert( ! TrustMe_ ) ;
  if ( TrustMe_ ) { 
    if ( CrsMatrixA_ == 0 ) AMESOS_CHK_ERR(10 );
    if( UseDataInPlace_ != 1 ) AMESOS_CHK_ERR( 10 ) ;
    if( Reindex_ )  AMESOS_CHK_ERR( 10 ) ;
    if( ! Problem_->GetLHS() )  AMESOS_CHK_ERR( 10 ) ;
    if( ! Problem_->GetRHS() )  AMESOS_CHK_ERR( 10 ) ;
    if( ! Problem_->GetLHS()->NumVectors() ) AMESOS_CHK_ERR( 10 ) ;
    if( ! Problem_->GetRHS()->NumVectors() ) AMESOS_CHK_ERR( 10 ) ; 
    SerialB_ = Problem_->GetRHS() ;
    SerialX_ = Problem_->GetLHS() ;
    NumVectors_ = SerialX_->NumVectors();
    if (MyPID_ == 0) {
      AMESOS_CHK_ERR(SerialX_->ExtractView(&SerialXBvalues_,&SerialXlda_ ));
      AMESOS_CHK_ERR(SerialB_->ExtractView(&SerialBvalues_,&SerialXlda_ ));
    }
  }
  
  
  PrivateParakleteData_->common_ = rcp(new paraklete_common());
  
  const Epetra_MpiComm* MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm());
  assert (MpiComm != 0);
  
  MPI_Comm PK_Comm;
  //
  //  Create an MPI group with MaxProcesses_ processes
  //
  if ( MaxProcesses_ != Comm().NumProc()) {
    if(ParakleteComm_)  {
      MPI_Comm_free(&ParakleteComm_);
      ParakleteComm_ = 0 ; 
    }
    std::vector<int> ProcsInGroup(MaxProcesses_);
    IamInGroup_ = false; 
    for (int i = 0 ; i < MaxProcesses_ ; ++i) {
      ProcsInGroup[i] = i;
      if ( Comm().MyPID() == i ) IamInGroup_ = true; 
    }
    
    MPI_Group OrigGroup, ParakleteGroup;
    MPI_Comm_group(MpiComm->GetMpiComm(), &OrigGroup);
    MPI_Group_incl(OrigGroup, MaxProcesses_, &ProcsInGroup[0], &ParakleteGroup);
    MPI_Comm_create(MpiComm->GetMpiComm(), ParakleteGroup, &ParakleteComm_);
    PK_Comm = ParakleteComm_ ; 
  } else {
    IamInGroup_ = true; 
    PK_Comm = MpiComm->GetMpiComm() ;
  }
  
  paraklete_common& pk_common =  *PrivateParakleteData_->common_ ;
  cholmod_common *cm = &(pk_common.cm) ;
  amesos_cholmod_l_start (cm) ;
  PK_DEBUG_INIT ("pk", cm) ;
  pk_common.nproc = MaxProcesses_ ;
  pk_common.myid = Comm().MyPID() ; 
  //pk_common.mpicomm = PK_Comm ; 
  cm->print = 1 ;
  cm->precise = TRUE ;
  cm->error_handler = my_handler ;
  
  pk_common.tol_diag = 0.001 ;
  pk_common.tol_offdiag = 0.1 ;
  pk_common.growth = 2. ;

  

  // "overhead" time for the following two methods is considered here
  AMESOS_CHK_ERR( ExportToSerial() );

  AMESOS_CHK_ERR( ConvertToParakleteCRS(true) );

  OverheadTime_ = AddTime("Total Amesos overhead time", OverheadTime_, 1);

  // All this time if PARAKLETE time
  AMESOS_CHK_ERR( PerformSymbolicFactorization() );

  NumSymbolicFact_++;

  IsSymbolicFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Paraklete::NumericFactorization() 
{
  if ( true || !TrustMe_ ) 
  { 
    IsNumericFactorizationOK_ = false;
    if (IsSymbolicFactorizationOK_ == false)
      AMESOS_CHK_ERR(SymbolicFactorization());
    
    ResetTimer(1); // "overhead" time

    Epetra_CrsMatrix *CrsMatrixA = dynamic_cast<Epetra_CrsMatrix *>(RowMatrixA_);
    if ( false && CrsMatrixA == 0 )   // hack to get around Bug #1502 
      AMESOS_CHK_ERR( CreateLocalMatrixAndExporters() ) ;
    assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );
    if ( AddZeroToDiag_ == 0 ) 
      assert( numentries_ == RowMatrixA_->NumGlobalNonzeros() );

    AMESOS_CHK_ERR( ExportToSerial() );
    
    if (  false && CrsMatrixA == 0 ) {  // continuation of hack to avoid bug #1502
      AMESOS_CHK_ERR( ConvertToParakleteCRS(true) );
    }  else {
      AMESOS_CHK_ERR( ConvertToParakleteCRS(false) );
    }

    OverheadTime_ = AddTime("Total Amesos overhead time", OverheadTime_, 1);
  }

  // this time is all for PARAKLETE
  AMESOS_CHK_ERR( PerformNumericFactorization() );

  NumNumericFact_++;

  IsNumericFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Paraklete::Solve() 
{
  Epetra_MultiVector* vecX=0;
  Epetra_MultiVector* vecB=0;

  if ( !TrustMe_ ) { 

    SerialB_ = Problem_->GetRHS() ;
    SerialX_ = Problem_->GetLHS() ;

    Epetra_MultiVector* OrigVecX ;
    Epetra_MultiVector* OrigVecB ;

    if (IsNumericFactorizationOK_ == false)
      AMESOS_CHK_ERR(NumericFactorization());
    
    ResetTimer(1);

    //
    //  Reindex the LHS and RHS 
    //
    OrigVecX = Problem_->GetLHS() ;
    OrigVecB = Problem_->GetRHS() ;
    
    if ( Reindex_ ) { 
#ifdef HAVE_AMESOS_EPETRAEXT
      vecX = StdIndexDomain_->StandardizeIndex( OrigVecX ) ;
      vecB = StdIndexRange_->StandardizeIndex( OrigVecB ) ;
#else
      AMESOS_CHK_ERR( -13 ) ; // Amesos_Paraklete can't handle non-standard indexing without EpetraExt 
#endif
    } else {
      vecX = OrigVecX ;
      vecB = OrigVecB ;
    } 
    
    if ((vecX == 0) || (vecB == 0))
      AMESOS_CHK_ERR(-1); // something wrong in input
    
    //  Extract Serial versions of X and B
    
    ResetTimer(0);

    //  Copy B to the serial version of B
    //
    if (false && UseDataInPlace_ == 1) {
      SerialB_ = vecB;
      SerialX_ = vecX;
      NumVectors_ = Problem_->GetRHS()->NumVectors() ; 
    } else {
      assert (UseDataInPlace_ == 0);
      
      if( NumVectors_ != Problem_->GetRHS()->NumVectors() ) {
	NumVectors_ = Problem_->GetRHS()->NumVectors() ; 
	SerialXextract_ = rcp( new Epetra_MultiVector(*SerialMap_,NumVectors_));
	SerialBextract_ = rcp (new Epetra_MultiVector(*SerialMap_,NumVectors_));
      }
      if (NumVectors_ != vecB->NumVectors())
	AMESOS_CHK_ERR(-1); // internal error 
      
      ImportRangeToSerial_ = rcp(new Epetra_Import ( *SerialMap_, vecB->Map() ) );
      if ( SerialBextract_->Import(*vecB,*ImportRangeToSerial_,Insert) )
	AMESOS_CHK_ERR( -1 ) ; // internal error
      
      SerialB_ = &*SerialBextract_ ;
      SerialX_ = &*SerialXextract_ ;
    }
    VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_, 0);

    //  Call PARAKLETE to perform the solve
    
    ResetTimer(0);
    if (MyPID_ == 0) {
      AMESOS_CHK_ERR(SerialB_->ExtractView(&SerialBvalues_,&SerialXlda_ ));
      AMESOS_CHK_ERR(SerialX_->ExtractView(&SerialXBvalues_,&SerialXlda_ ));
      if (SerialXlda_ != NumGlobalElements_)
	AMESOS_CHK_ERR(-1);
    }

    OverheadTime_ = AddTime("Total Amesos overhead time", OverheadTime_, 1);
  }
  if ( MyPID_ == 0) {
    if ( NumVectors_ == 1 ) {
      for ( int i = 0 ; i < NumGlobalElements_ ; i++ ) 
	SerialXBvalues_[i] = SerialBvalues_[i] ;
    } else {
      SerialX_->Scale(1.0, *SerialB_ ) ;    // X = B (Klu overwrites B with X)
    }
  }
  if ( IamInGroup_ ) 
    for (int i = 0; i < NumVectors_ ; i++ ) { 
      amesos_paraklete_solve(  &*PrivateParakleteData_->LUnumeric_, &*PrivateParakleteData_->LUsymbolic_,
			&SerialXBvalues_[i*SerialXlda_],  &*PrivateParakleteData_->common_ );
    }
  SolveTime_ = AddTime("Total solve time", SolveTime_, 0);

  //  Copy X back to the original vector

  ResetTimer(0);
  ResetTimer(1);

  if (UseDataInPlace_ == 0) {
    ImportDomainToSerial_ = rcp(new Epetra_Import ( *SerialMap_, vecX->Map() ) );
    vecX->Export( *SerialX_, *ImportDomainToSerial_, Insert ) ;

  } // otherwise we are already in place.

  VecRedistTime_ = AddTime("Total vector redistribution time", VecRedistTime_, 0);

#if 0
  //
  //  ComputeTrueResidual causes TestOptions to fail on my linux box //  Bug #1417
  if (ComputeTrueResidual_)
    ComputeTrueResidual(*SerialMatrix_, *vecX, *vecB, UseTranspose(), "Amesos_Paraklete");
#endif

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Paraklete");

  OverheadTime_ = AddTime("Total Amesos overhead time", OverheadTime_, 1);

  ++NumSolve_;

  return(0) ;
}

// ================================================ ====== ==== ==== == =

void Amesos_Paraklete::PrintStatus() const
{

  if (MyPID_) return;

  PrintLine();

  std::cout << "Amesos_Paraklete : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << numentries_ << " nonzeros" << std::endl;
  std::cout << "Amesos_Paraklete : Nonzero elements per row = "
       << 1.0*numentries_/NumGlobalElements_ << std::endl;
  std::cout << "Amesos_Paraklete : Percentage of nonzero elements = "
       << 100.0*numentries_/(pow(double(NumGlobalElements_),double(2.0))) << std::endl;
  std::cout << "Amesos_Paraklete : Use transpose = " << UseTranspose_ << std::endl;

  PrintLine();

  return;

}

// ================================================ ====== ==== ==== == =

void Amesos_Paraklete::PrintTiming() const
{
  if (MyPID_) return;

  double ConTime = GetTime(MtxConvTime_);
  double MatTime = GetTime(MtxRedistTime_);
  double VecTime = GetTime(VecRedistTime_);
  double SymTime = GetTime(SymFactTime_);
  double NumTime = GetTime(NumFactTime_);
  double SolTime = GetTime(SolveTime_);
  double OveTime = GetTime(OverheadTime_);

  if (NumSymbolicFact_)
    SymTime /= NumSymbolicFact_;

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  std::string p = "Amesos_Paraklete : ";
  PrintLine();

  std::cout << p << "Time to convert matrix to Paraklete format = "
       << ConTime << " (s)" << std::endl;
  std::cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << std::endl;
  std::cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << std::endl;
  std::cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << std::endl;
  std::cout << p << "Time for sym fact = "
       << SymTime * NumSymbolicFact_ << " (s), avg = " << SymTime << " (s)" << std::endl;
  std::cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << std::endl;
  std::cout << p << "Time for num fact = "
       << NumTime * NumNumericFact_ << " (s), avg = " << NumTime << " (s)" << std::endl;
  std::cout << p << "Number of solve phases = "
       << NumSolve_ << std::endl;
  std::cout << p << "Time for solve = "
       << SolTime * NumSolve_ << " (s), avg = " << SolTime << " (s)" << std::endl;

  double tt = SymTime * NumSymbolicFact_ + NumTime * NumNumericFact_ + SolTime * NumSolve_;
  if (tt != 0)
  {
    std::cout << p << "Total time spent in Amesos = " << tt << " (s) " << std::endl;
    std::cout << p << "Total time spent in the Amesos interface = " << OveTime << " (s)" << std::endl;
    std::cout << p << "(the above time does not include PARAKLETE time)" << std::endl;
    std::cout << p << "Amesos interface time / total time = " << OveTime / tt << std::endl;
  }

  PrintLine();

  return;
}


