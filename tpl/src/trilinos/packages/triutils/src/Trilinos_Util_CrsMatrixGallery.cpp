// @HEADER
// ***********************************************************************
//
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Trilinos_Util.h"
#include <string>

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace std;

const double UNDEF = -99999.87;
const bool Scaling = false;

// ================================================ ====== ==== ==== == =
Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery(const string name,
    const Epetra_Comm & comm,
    bool UseLongLong) :
  comm_(&comm), name_(name), UseLongLong_(UseLongLong)
{
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(!UseLongLong)
    throw "Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery: UseLongLong = false but no 32 bit global indices";
#endif
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong)
    throw "Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery: UseLongLong = true but no 64 bit global indices";
#endif

  ZeroOutData();
  // verbosity level (now always false)
  // if( comm_->MyPID()==0 ) verbose_ = true;
  verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR [CrsMatrixGallery]: ";
  OutputMsg = "CrsMatrixGallery: ";

}

// ================================================ ====== ==== ==== == =
Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery(const string name,
    const Epetra_Map & map ) :
  comm_(&(map.Comm())), name_(name)
{
  ZeroOutData();
  // verbosity level
  if( comm_->MyPID()==0 ) verbose_ = true;
  else verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR [Trilinos_Util::CrsMatrixGallery]: ";
  OutputMsg = "Trilinos_Util::CrsMatrixGallery: ";

  map_ = new Epetra_Map(map);
  UseLongLong_ = map.GlobalIndicesLongLong();
  NumGlobalElements_ = map_->NumGlobalElements64();
  NumMyElements_ = map_->NumMyElements();
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(map_->GlobalIndicesInt()) {
    MyGlobalElements_int_ = map_->MyGlobalElements();
    UseLongLong_ = false;
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    if(map_->GlobalIndicesLongLong()) {
      MyGlobalElements_LL_ = map_->MyGlobalElements64();
      UseLongLong_ = true;
    }
    else
#endif
      throw "Trilinos_Util::CrsMatrixGallery::CrsMatrixGallery: Global Indices unknown";
}

// ================================================ ====== ==== ==== == =
Trilinos_Util::CrsMatrixGallery::~CrsMatrixGallery(void)
{

  // linear problem
  if( LinearProblem_ != NULL ) delete LinearProblem_;

  // Crs data
  if( matrix_ != NULL ) delete matrix_;
  if( ExactSolution_ != NULL ) delete ExactSolution_;
  if( StartingSolution_ != NULL ) delete StartingSolution_;
  if( rhs_ != NULL ) delete rhs_;
  if( map_ != NULL ) delete map_;

  // vectors
  if( VectorA_ != NULL ) delete VectorA_;
  if( VectorB_ != NULL ) delete VectorB_;
  if( VectorC_ != NULL ) delete VectorC_;
  if( VectorD_ != NULL ) delete VectorD_;
  if( VectorE_ != NULL ) delete VectorE_;
  if( VectorF_ != NULL ) delete VectorF_;
  if( VectorG_ != NULL ) delete VectorG_;

  // put to default values
  ZeroOutData();
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const int value)
{

  if( parameter == "problem_size" ) {
    if( value <= 0 ) {
      cerr << ErrorMsg << "problem size must be greater than 1\n";
      return -1;
    }
    if( map_ != NULL ) {
      cerr << ErrorMsg << "map object already set. Continuing with\n"
        << ErrorMsg << "problemSize = " << NumGlobalElements_ << endl;
      return -2;
    }
    NumGlobalElements_ = value;
    return 0;

  } else if( parameter == "nx" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "nx must be greater than 0\n";
      return -1;
    }

    nx_ = value;
    return 0;

  } else if( parameter == "ny" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "ny must be greater than 0\n";
      return -1;
    }

    ny_ = value;
    return 0;

  } else if( parameter == "nz" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "nz must be greater than 0\n";
      return -1;
    }

    nz_ = value;
    return 0;
  } else if( parameter == "mx" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "mx must be greater than 0\n";
      return -1;
    }

    mx_ = value;
    return 0;

  } else if( parameter == "my" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "my must be greater than 0\n";
      return -1;
    }

    my_ = value;
    return 0;

  } else if( parameter == "mz" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "mz must be greater than 0\n";
      return -1;
    }

    mz_ = value;
    return 0;

  } else if( parameter == "num_pde_eqns" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "num pde eqns must be greater than 0\n";
      return -1;
    }

    NumPDEEqns_ = value;
    return 0;

  } else if( parameter == "num_vectors" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "num_vectors must be greater than 0\n";
      return -1;
    }

    NumVectors_ = value;
    return 0;
  } else if( parameter == "output" ) {

    if( value != 0 && value != 1) {
      cerr << ErrorMsg << "output level should be 0 or 1" << endl;
      return -1;
    }

    if (value == 0)
      verbose_ = false;
    if (value == 1) {
      if (comm_->MyPID() == 0)
        verbose_ = true;
    }

    return 0;
  }

  cerr << ErrorMsg << "input string (" << parameter << ") not valid\n";
  return -2;

}

// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const string value )
{

  if( parameter == "problem_type" ) {
    name_ = value;
  }
  else if( parameter == "map_type" ) {
    MapType_ = value;
  }
  else if( parameter == "exact_solution" ) {
    ExactSolutionType_ = value;
  }
  else if( parameter == "matrix_name" ) {
    FileName_ = value;
  }
  else if( parameter == "starting_solution" ) {
    StartingSolutionType_ = value;
  }
  else if( parameter == "rhs_type" ) {
    RhsType_ = value;
  }
  else if( parameter == "noncontiguos_map" ) {
    ContiguousMap_ = false;
  }
  else if( parameter == "output" ) {
    if( value == "none" ) verbose_ = false;
    else {
      if( value == "proc 0" ) {
        if( comm_->MyPID()==0 ) verbose_ = true;
        else verbose_ = false;
      } else {
        verbose_ = true;
      }
    }
  } else if( parameter == "expand_type" ) {
    ExpandType_ = value;
  } else {
    cerr << ErrorMsg << "wrong input parameter (" << parameter << ")\n";
    return -1;
  }

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const double value)
{

  if( parameter == "a" ) {
    a_ = value;
    return 0;
  } else if( parameter == "b" ) {
    b_ = value;
    return 0;
  } else if( parameter == "c" ) {
    c_ = value;
    return 0;
  } else if( parameter == "d" ) {
    d_ = value;
    return 0;
  } else if( parameter == "e" ) {
    e_ = value;
    return 0;
  } else if( parameter == "f" ) {
    f_ = value;
    return 0;
  } else if( parameter == "g" ) {
    g_ = value;
    return 0;
  } else if( parameter == "conv" ) {
    conv_ = value;
    return 0;
  } else if( parameter == "diff" ) {
    diff_ = value;
    return 0;
  } else if( parameter == "source" ) {
    source_ = value;
    return 0;
  } else if( parameter == "alpha" ) {
    alpha_ = value;
    return 0;
  } else if( parameter == "epsilon" ) {
    epsilon_= value;
    return 0;
  } else if( parameter == "lx" ) {
    lx_ = value;
    return 0;
  } else if( parameter == "ly" ) {
    ly_ = value;
    return 0;
  } else if( parameter == "lz" ) {
    lz_ = value;
    return 0;
  }

  cerr << ErrorMsg << "input string not valid\n";
  return -2;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::Set(const string parameter, const Epetra_Vector & value)
{

  if( value.Map().SameAs(*map_) == false ) {
    cerr << ErrorMsg << "input vector must have the same map used to\n"
      << ErrorMsg << "create the Trilinos_Util::CrsMatrixGallery object. Continuing\n";
    return -2;
  }

  if( parameter == "a" ) {
    VectorA_ = new Epetra_Vector(value);
  } else if( parameter == "b" ) {
    VectorB_ = new Epetra_Vector(value);
  }
  else if( parameter == "c" ) {
    VectorC_ = new Epetra_Vector(value);
  }
  else if( parameter == "d" ) {
    VectorD_ = new Epetra_Vector(value);
  }
  else if( parameter == "e" ) {
    VectorE_ = new Epetra_Vector(value);
  }
  else if( parameter == "f" ) {
    VectorF_ = new Epetra_Vector(value);
  }
  else if( parameter == "g" ) {
    VectorG_ = new Epetra_Vector(value);
  } else {
    cerr << ErrorMsg << "input string not valid\n";
    return -3;
  }

  return 0;
}

// ================================================ ====== ==== ==== == =

int Trilinos_Util::CrsMatrixGallery::Set(Trilinos_Util::CommandLineParser & CLP)
{
  int count;

  string Options[15];

  // all options with strings
  count = 0;
  Options[count++] = "problem_type";
  Options[count++] = "map_type";
  Options[count++] = "exact_solution";
  Options[count++] = "matrix_name";
  Options[count++] = "starting_solution";
  Options[count++] = "output";
  Options[count++] = "expand_type";
  Options[count++] = "rhs_type";

  for( int i=0 ; i<count ; i++ ) {
    string parameter = "-"+Options[i];
    if( CLP.Has(parameter) ) {
      string value = CLP.Get(parameter,"not-set");
      Set(Options[i],value);

    }
  }

  // all options with integers
  Options[0] = "problem_size";
  Options[1] = "nx";
  Options[2] = "ny";
  Options[3] = "nz";
  Options[4] = "mx";
  Options[5] = "my";
  Options[6] = "mz";
  Options[7] = "num_pde_eqns";

  for(  int i=0 ; i<8 ; i++ ) {
    string parameter = "-"+Options[i];
    if( CLP.Has(parameter) ) {
      Set(Options[i],CLP.Get(parameter,(int)1));
    }
  }

  // all options with doubles
  Options[0]  = "a";
  Options[1]  = "b";
  Options[2]  = "c";
  Options[3]  = "d";
  Options[4]  = "e";
  Options[5]  = "f";
  Options[6]  = "g";
  Options[7]  = "conv";
  Options[8]  = "diff";
  Options[9]  = "source";
  Options[10] = "alpha";
  Options[11] = "lx";
  Options[12] = "ly";
  Options[13] = "lz";

  for( int i=0 ; i<14 ; i++ ) {
    string parameter = "-"+Options[i];

    if( CLP.Has(parameter) ) {
      Set(Options[i],CLP.Get(parameter,1.0));
    }

  }

  return 0;
}

// ================================================ ====== ==== ==== == =
Epetra_CrsMatrix * Trilinos_Util::CrsMatrixGallery::GetMatrix(void)
{
  if( matrix_ == NULL ) CreateMatrix();
  return( matrix_ );
}

// ================================================ ====== ==== ==== == =
Epetra_CrsMatrix & Trilinos_Util::CrsMatrixGallery::GetMatrixRef(void)
{
  if( matrix_ == NULL ) CreateMatrix();
  return( *matrix_ );
}

// ================================================ ====== ==== ==== == =
Epetra_MultiVector * Trilinos_Util::CrsMatrixGallery::GetExactSolution(void)
{
  if( ExactSolution_ == NULL ) CreateExactSolution();
  return ExactSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_MultiVector * Trilinos_Util::CrsMatrixGallery::GetStartingSolution(void)
{
  if( StartingSolution_ == NULL ) CreateStartingSolution();
  return StartingSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_MultiVector * Trilinos_Util::CrsMatrixGallery::GetRHS(void)
{
  if( rhs_ == NULL ) CreateRHS();
  return rhs_;
}

// ================================================ ====== ==== ==== == =
const Epetra_Map * Trilinos_Util::CrsMatrixGallery::GetMap(void)
{
  if( map_ == NULL ) CreateMap();

  return map_;
}

// ================================================ ====== ==== ==== == =
const Epetra_Map & Trilinos_Util::CrsMatrixGallery::GetMapRef(void)
{
  if( map_ == NULL ) CreateMap();

  return *map_;
}

// ================================================ ====== ==== ==== == =
Epetra_LinearProblem * Trilinos_Util::CrsMatrixGallery::GetLinearProblem(void)
{
  // pointers, not really needed
  Epetra_CrsMatrix * A;
  Epetra_MultiVector * RHS;
  Epetra_MultiVector * StartingSolution;

  A = GetMatrix();
  RHS = GetRHS();
  StartingSolution = GetStartingSolution();

  // create linear problem
  if( LinearProblem_ != NULL ) delete LinearProblem_;
  LinearProblem_ = new Epetra_LinearProblem(A,StartingSolution,RHS);

  return LinearProblem_;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ComputeResidual(double* residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution_ are
  //  created by CreateRHS if necessary)
  if( rhs_ == NULL ) CreateRHS();

  Epetra_MultiVector Ax(*map_, NumVectors_);

  matrix_->Multiply(false, *StartingSolution_, Ax);
  Ax.Update(1.0, *rhs_, -1.0);
  Ax.Norm2(residual);

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ComputeDiffBetweenStartingAndExactSolutions(double* residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( rhs_ == NULL ) CreateRHS();

  Epetra_MultiVector temp(*map_, NumVectors_);

  temp.Update(1.0, *ExactSolution_, -1.0, *StartingSolution_, 0.0);
  temp.Norm2(residual);

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::TCreateMap(void)
{

  Epetra_Time Time(*comm_);

  if( verbose_ ) {
    cout << OutputMsg << "Creating Map `" << MapType_ << "'...\n";
  }

  // first get the problem size. For some problems. the user can
  // specify the problem size using different parameters (e.g.,
  // nx and ny for a 2D Laplace problem). I need the internal
  // variable NumGlobalElements_ properly set before continuing.
  // NOTE: for HB problems, this value has already been set

  if( name_ == "diag" || name_ == "tridiag"  ||
      name_ == "laplace_1d" || name_ == "laplace_1d_n" ||
      name_ == "eye" ||
      name_ == "lehmer" || name_ == "minij" ||
      name_ == "ris" || name_ == "hilbert" ||
      name_ == "jordblock" || name_ == "cauchy" ||
      name_ == "fielder" || name_ == "hanowa" ||
      name_ == "kms" || name_ == "parter" ||
      name_ == "pei" || name_ == "ones" ||
      name_ == "vander" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 ) NumGlobalElements_ = nx_;
      else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            true,
            ErrorMsg << "problem size not correct (" << NumGlobalElements_ << ")"
            );
      }
    }
  }

  else if( name_ == "laplace_2d" || name_ == "laplace_2d_n"
      || name_ == "laplace_2d_bc"
      || name_ == "cross_stencil_2d"
      || name_ == "laplace_2d_9pt" || name_ == "recirc_2d"
      || name_ == "uni_flow_2d" || name_ == "recirc_2d_divfree"
      || name_ == "stretched_2d" ) {

    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 && ny_ > 0 )
        NumGlobalElements_ = nx_*ny_;
      else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            true,
            ErrorMsg << "Problem size not correct (" << NumGlobalElements_ << ")\n"
            << ErrorMsg << "It should be a perfect square"
            );
      }
    }

    if( verbose_ ) {
      cout << OutputMsg << "nx = " << nx_ << ", ny = " << ny_ << endl;
    }

  }

  else if( name_ == "laplace_3d" || name_ == "cross_stencil_3d" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 && ny_ > 0 && nz_ > 0 )
        NumGlobalElements_ = nx_*ny_*nz_;
      else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            true,
            ErrorMsg << "Problem size not correct (" << NumGlobalElements_  << ")\n"
            << ErrorMsg << "It should be a perfect cube"
            );
      }
    }

    if( verbose_ ) {
      cout << OutputMsg << "nx = " << nx_ << ", ny = " << ny_ << ", nz = " << nz_ << endl;
    }

  } else if( name_ == "hb" || name_ == "matrix_market" ||
      name_ == "triples_sym" || name_ == "triples_nonsym" ) {
    // The global number of elements has been set in reading matrix
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        NumGlobalElements_ <= 0,
        ErrorMsg << "Problem size not correct (" << NumGlobalElements_ << ")"
        );

  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        ErrorMsg << "matrix name is incorrect or not set ("
        << name_ << ")"
        );
  }

  std::vector<int_type>& MapMap = MapMapRef<int_type>();

  if (! ContiguousMap_ ) {
    //
    //  Populate MapMap[] with NumGlobalElements_ numbers randomly
    //  chosen from the set of integers ranging from 0 to 2*NumGlobalElements-1
    //
    MapMap.resize( NumGlobalElements_ ) ;
    Epetra_IntSerialDenseVector sortable_values(2*NumGlobalElements_) ;
    sortable_values.Random();
    Epetra_IntSerialDenseVector sortable_positions(2*NumGlobalElements_) ;
    for( int_type i =0 ; i < 2*NumGlobalElements_; i++ ) {
      sortable_positions[i] = i;
    }
    Epetra_Util Utils;
    int *ArrayOfIntPointers[1];
    //    double *ArrayOfDoublePointers[1];
    ArrayOfIntPointers[0] =  &sortable_positions[0];
    //    ArrayOfDoublePointers[0] =  &sortable_values[0];
    Utils.Sort( true, NumGlobalElements_*2, &sortable_values[0],
        0, (double**)0, 1, &ArrayOfIntPointers[0] );

    for( int_type i =0 ; i < NumGlobalElements_; i++ ) {
      MapMap[i] = sortable_positions[i];
    }
    //
    //  Make sure that all processes have the same indices in MapMap
    //
    comm_->Broadcast( &MapMap[0], NumGlobalElements_, 0 ) ;
  }
  // check out whether one is using only one proc or not.
  // If yes, creation of map is straightforward. Then return.

  if( comm_->NumProc() == 1 ) {

    if (ContiguousMap_ )
      map_ = new Epetra_Map((int_type) NumGlobalElements_,(int_type) 0,*comm_);
    else
      map_ = new Epetra_Map((int_type) NumGlobalElements_,NumMyElements_,&MapMap[0], (int_type) 0,*comm_);

  } else {

    // Here below more than one processor.

    if( MapType_ == "linear" ) {

      map_ = new Epetra_Map ((int_type) NumGlobalElements_,(int_type) 0,*comm_);
      if ( ! ContiguousMap_ ) {
        //
        //  map_ gives us NumMyElements and MyFirstElement;
        //
        int NumMyElements = map_->NumMyElements();
        int_type MyFirstElement = map_->MinMyGID64();
        assert( MyFirstElement+NumMyElements ==  (int_type) map_->MaxMyGID64());
        delete map_;
        map_ = new Epetra_Map( NumGlobalElements_, NumMyElements, &MapMap[MyFirstElement], 0, *comm_);
      }

    } else if( MapType_ == "box" ) {

      if( mx_ == -1 || my_ == -1 ) {
        mx_ = (int)sqrt((double)(comm_->NumProc()));
        my_ = mx_;

        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            mx_ * my_ != comm_->NumProc(),
            ErrorMsg << "number of processes must be perfect square\n"
            << ErrorMsg << "otherwise set mx and my\n"
            );
      } else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            mx_ * my_ != comm_->NumProc(),
            ErrorMsg << "mx*my != number of processes ("
            << mx_ * my_ << " != " << comm_->NumProc()  << ")"
            );
      }

      if( verbose_ ) {
        cout << OutputMsg << "mx = " << mx_ << ", my = " << my_ << endl;
      }

      SetupCartesianGrid2D();

      // how to divide the axis

      int modx = (nx_+(nx_%mx_))/mx_;
      int mody = (ny_+(ny_%my_))/my_;

      int MyPID = comm_->MyPID(), startx, starty, endx, endy;
      int xpid = MyPID%mx_;
      int ypid = MyPID/mx_;

      startx = xpid*modx;
      if( (xpid+1)*modx < nx_ ) endx = (xpid+1)*modx;
      else endx = nx_;
      starty = ypid*mody;
      if( (ypid+1)*mody < ny_ ) endy = (ypid+1)*mody;
      else endy = ny_;

      int NumMyElements = (endx-startx)*(endy-starty);
      int_type * MyGlobalElements = new int_type[NumMyElements];
      int count = 0;

      for( int i=startx ; i<endx ; ++i ) {
        for( int j=starty ; j<endy ; ++j ) {
          MyGlobalElements[count++] = i+((int_type)j)*nx_;
        }
      }

      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);

      // I delete this guy so that this case is not different from the
      // others, and I don't have to clean up this mess while
      // destroying the object.

      delete [] MyGlobalElements;

    } else if( MapType_ == "cube" ) {

      if( mx_ == -1 || my_ == -1 || mz_ == -1 ) {
        mx_ = (int)pow((double)(comm_->NumProc()),0.333334);
        my_ = mx_;
        mz_ = mx_;

        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            mx_ * my_ * mz_ != comm_->NumProc(),
            ErrorMsg << "number of processes must be perfect cube\n"
            << ErrorMsg << "otherwise set mx, my, and mz"
            );
      } else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            mx_ * my_ * mz_ != comm_->NumProc(),
            ErrorMsg << "mx*my*mz != number of processes ("
            << mx_ * my_ * mz_ << " != " << comm_->NumProc()
            << ")"
            );
      }

      if( verbose_ ) {
        cout << OutputMsg << "mx = " << mx_ << ", my = " << my_ << ", mz = " << mz_ << endl;
      }

      SetupCartesianGrid3D();

      // how to divide the axis

      int modx = (nx_+(nx_%mx_))/mx_;
      int mody = (ny_+(ny_%my_))/my_;
      int modz = (nz_+(nz_%mz_))/mz_;

      int MyPID = comm_->MyPID(), startx, starty, startz, endx, endy, endz;
      int mxy  = mx_*my_;
      int zpid = MyPID/mxy;
      int xpid = (MyPID%mxy)%mx_;
      int ypid = (MyPID%mxy)/mx_;

      startx = xpid*modx;
      if( (xpid+1)*modx < nx_ ) endx = (xpid+1)*modx;
      else endx = nx_;
      starty = ypid*mody;
      if( (ypid+1)*mody < ny_ ) endy = (ypid+1)*mody;
      else endy = ny_;
      startz = zpid*modz;
      if( (zpid+1)*modz < nz_ ) endz = (zpid+1)*modz;
      else endz = nz_;

      int NumMyElements = (endx-startx)*(endy-starty)*(endz-startz);
      int_type * MyGlobalElements = new int_type[NumMyElements];
      int count = 0;

      for( int i=startx ; i<endx ; ++i ) {
        for( int j=starty ; j<endy ; ++j ) {
          for( int k=startz ; k<endz ; ++k ) {
            MyGlobalElements[count++] = i+((int_type)j)*nx_+((int_type)k)*(((int_type)nx_)*ny_);
          }
        }
      }

      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);

      // I delete this guy so that this case is not different from the
      // others, and I don't have to clean up this mess while
      // destroying the object.

      delete [] MyGlobalElements;

    } else if( MapType_ == "interlaced" ) {

      // this is the first funky map. Nodes are assigned so that
      // node 0 is given to proc 0, node 1 to proc 1, and
      // node i to proc i%NumProcs. Probably not the best, but it
      // results in decompositions with lots of boundary nodes.

      int NumProcs = comm_->NumProc();
      int MyPID = comm_->MyPID();

      int NumMyElements = NumGlobalElements_/NumProcs;
      if( MyPID < NumGlobalElements_%NumProcs ) NumMyElements++;

      int count = 0;
      int_type * MyGlobalElements = new int_type[NumMyElements];

      for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) {
        if( i%NumProcs == MyPID )
          MyGlobalElements[count++] = i;
      }

      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          count != NumMyElements,
          ErrorMsg << "something went wrong in CreateMap\n"
          << ErrorMsg << "count = " << count << ", NumMyElements = "
          << NumMyElements
          );

      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);
      delete [] MyGlobalElements;

    } else if( MapType_ == "random" ) {

      // this is even funkier. Random decomposition of nodes into procs.
      // It should result in a more ordered decomposition than "interlaced"
      // This is the idea: I create the map on proc 0, then I broadcast
      // it to all procs. This is not very efficient, but saves some
      // MPI calls.

      int * part = new int[NumGlobalElements_];

      if( comm_->MyPID() == 0 ) {
        Epetra_Util Util;

        for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) {
          unsigned int r = Util.RandomInt();
          part[i] = r%(comm_->NumProc());
        }
      }

      comm_->Broadcast(part,NumGlobalElements_,0);

      // count the elements assigned to this proc
      int NumMyElements = 0;
      for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) {
        if( part[i] == comm_->MyPID() ) NumMyElements++;
      }

      // get the loc2global list
      int_type * MyGlobalElements = new int_type[NumMyElements];
      int count = 0;
      for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) {
        if( part[i] == comm_->MyPID() ) MyGlobalElements[count++] = i;
      }

      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,
          0,*comm_);

      delete [] MyGlobalElements;
      delete [] part;

    } else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg << "MapType has an incorrect value (" << MapType_ << ")"
          );
    }
  }

  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  map_->MyGlobalElementsPtr(MyGlobalElementsPtr<int_type>());

  if( verbose_ ) {
    cout << OutputMsg << "Time to create Map: "
      << Time.ElapsedTime() << " (s)\n";
  }

  return;
}

void Trilinos_Util::CrsMatrixGallery::CreateMap(void)
{
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_)
    TCreateMap<long long>();
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_)
      TCreateMap<int>();
    else
#endif
      throw "Trilinos_Util::CrsMatrixGallery::CreateMap: failed";
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::TCreateMatrix(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating Matrix...\n";
  }

  // HB matrices are different, as their dimension has to be read before.
  // Here the idea is to read the matrix on proc 0, then build the
  // map, then redistribute it linearly

  if( name_ == "hb" || name_ == "matrix_market" ||
      name_ == "triples_sym" || name_ == "triples_nonsym" ) {

    Epetra_Time Time(*comm_);
    TReadMatrix<int_type>();
    if( verbose_ ) {
      cout << OutputMsg << "Time to create matrix: "
        << Time.ElapsedTime() << " (s)\n";
    }

  } else {

    if( map_ == NULL ) CreateMap();

    Epetra_Time Time(*comm_);

    if( name_ == "diag" ) CreateMatrixDiag<int_type>();

    else if( name_ == "eye" ) CreateEye<int_type>();

    else if( name_ == "tridiag" ) CreateMatrixTriDiag<int_type>();

    else if( name_ == "laplace_1d" ) CreateMatrixLaplace1d<int_type>();

    else if( name_ == "laplace_1d_n" ) CreateMatrixLaplace1dNeumann<int_type>();

    else if( name_ == "laplace_2d" ) CreateMatrixLaplace2d<int_type>();

    else if( name_ == "laplace_2d_bc" ) CreateMatrixLaplace2d_BC<int_type>();

    else if( name_ == "laplace_2d_n" ) CreateMatrixLaplace2dNeumann<int_type>();

    else if( name_ == "laplace_2d_9pt" ) CreateMatrixLaplace2d_9pt<int_type>();

    else if( name_ == "stretched_2d" ) CreateMatrixStretched2d<int_type>();

    else if( name_ == "recirc_2d" ) CreateMatrixRecirc2d<int_type>();

    else if( name_ == "recirc_2d_divfree" ) CreateMatrixRecirc2dDivFree<int_type>();

    else if( name_ == "uni_flow_2d" ) CreateMatrixUniFlow2d<int_type>();

    else if( name_ == "laplace_3d" ) CreateMatrixLaplace3d<int_type>();

    else if( name_ == "cross_stencil_2d" ) CreateMatrixCrossStencil2d<int_type>();

    else if( name_ == "cross_stencil_3d" ) CreateMatrixCrossStencil3d<int_type>();

    else if( name_ == "lehmer" ) CreateMatrixLehmer<int_type>();

    else if( name_ == "minij" ) CreateMatrixMinij<int_type>();

    else if( name_ == "ris" ) CreateMatrixRis<int_type>();

    else if( name_ == "hilbert" ) CreateMatrixHilbert<int_type>();

    else if( name_ == "jordblock" ) CreateMatrixJordblock<int_type>();

    else if( name_ == "cauchy" ) CreateMatrixCauchy<int_type>();

    else if( name_ == "fiedler" ) CreateMatrixFiedler<int_type>();

    else if( name_ == "hanowa" ) CreateMatrixHanowa<int_type>();

    else if( name_ == "kms" ) CreateMatrixKMS<int_type>();

    else if( name_ == "parter" ) CreateMatrixParter<int_type>();

    else if( name_ == "pei" ) CreateMatrixPei<int_type>();

    else if( name_ == "ones" ) CreateMatrixOnes<int_type>();

    else if( name_ == "vander" ) CreateMatrixVander<int_type>();

    else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg << "matrix name is incorrect or not set ("
          << name_ << ")"
          );
    }

    if( verbose_ ) {
      cout << OutputMsg << "Time to create matrix: "
        << Time.ElapsedTime() << " (s)\n";
    }
  }

  matrix_->OptimizeStorage();

  return;
}

void Trilinos_Util::CrsMatrixGallery::CreateMatrix(void)
{
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_)
    TCreateMatrix<long long>();
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_)
      TCreateMatrix<int>();
    else
#endif
      throw "Trilinos_Util::CrsMatrixGallery::CreateMatrix: failed";
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::TCreateExactSolution(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating exact solution `"
      << ExactSolutionType_ << "'...\n";
  }

  if( map_ == NULL ) CreateMap();

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  if( ExactSolution_ == NULL ) {

    ExactSolution_ = new Epetra_MultiVector(*map_, NumVectors_);

    if( ExactSolutionType_ == "random" ) {

      ExactSolution_->Random();

    } else if( ExactSolutionType_ == "constant" ) {

      ExactSolution_->PutScalar(1.0);

    } else if( ExactSolutionType_ == "quad_x" ) {

      int_type*& MyGlobalElements2 = MyGlobalElementsPtr<int_type>();

      // always suppose to have Dirichlet boundary
      // conditions, and those points have already been eliminated
      // from the matrix
      double hx = lx_/(NumGlobalElements_+1);
      for( int i=0 ; i<NumMyElements_ ; i++ ) {
        double x = (MyGlobalElements2[i]+1)*hx;
        for (int j = 0 ; j < NumVectors_ ; ++j)
          (*ExactSolution_)[j][i] = x*(1.-x);
      }

    } else if( ExactSolutionType_ == "quad_xy" ) {

      SetupCartesianGrid2D();

      double hx = lx_/(nx_+1);
      double hy = ly_/(ny_+1);

      for( int i=0 ; i<NumMyElements_ ; ++i ) {
        int ix, iy;
        ix = (MyGlobalElements[i])%nx_;
        iy = (MyGlobalElements[i] - ix)/nx_;
        double x = hx*(ix+1);
        double y = hy*(iy+1);
        double u;
        ExactSolQuadXY(x,y,u);

        for (int j = 0 ; j < NumVectors_ ; ++j)
          (*ExactSolution_)[j][i] = u;

      }


    } else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg << "exact solution type is not correct : "
          << ExactSolutionType_ << "\n"
          << ErrorMsg << "It should be:\n"
          << ErrorMsg << "<random> / <constant> / <quad_x> / <quad_xy>"
          );
    }
  }

  return;

}

void Trilinos_Util::CrsMatrixGallery::CreateExactSolution(void)
{
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_)
    TCreateExactSolution<long long>();
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_)
      TCreateExactSolution<int>();
    else
#endif
      throw "Trilinos_Util::CrsMatrixGallery::CreateExactSolution: failed";
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::CreateStartingSolution(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating starting solution `"
      << StartingSolutionType_ << "'...\n";
  }

  if( map_ == NULL ) CreateMap();

  if( StartingSolution_ == NULL ) {
    StartingSolution_ = new Epetra_MultiVector(*map_, NumVectors_);
    if( StartingSolutionType_ == "random" ) {
      StartingSolution_->Random();
    } else if( StartingSolutionType_ == "zero" ) {
      StartingSolution_->PutScalar(0.0);
    } else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg << "starting solution type is not correct : "
          << StartingSolutionType_
          );
    }
  }

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::TCreateRHS(void)
{

  if( map_ == NULL ) CreateMap();
  if( matrix_ == NULL ) CreateMatrix();
  if( ExactSolution_ == NULL )  CreateExactSolution();

  if( rhs_ != NULL ) delete rhs_;

  Epetra_Time Time(*comm_);

  if( verbose_ ) {
    cout << OutputMsg << "Creating RHS `" << RhsType_ << "' ...\n";
  }

  rhs_ = new Epetra_MultiVector(*map_,NumVectors_);
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  if( RhsType_ == "from_exact_solution" ) {

    matrix_->Multiply(false,*ExactSolution_,*rhs_);

  } else if( RhsType_ == "exact_rhs_uni_flow_2d" ) {

    // need to set a_ and b_ too
    if( conv_ == UNDEF ) conv_ = 1;
    if( diff_ == UNDEF ) diff_ = 1e-5;
    if( alpha_ == UNDEF ) alpha_ = 1e-5;

    SetupCartesianGrid2D();

    double hx = lx_/(nx_+1);
    double hy = ly_/(ny_+1);

    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int ix, iy;
      ix = (MyGlobalElements[i])%nx_;
      iy = (MyGlobalElements[i] - ix)/nx_;
      double x = hx*(ix+1);
      double y = hy*(iy+1);
      double u, ux, uy, uxx, uyy;
      ExactSolQuadXY(x,y,u,ux,uy,uxx,uyy);

      for (int j = 0 ; j < NumVectors_ ; ++j) {
        (*rhs_)[j][i] = -diff_*( uxx + uyy )  // -b_ \nabla u
          + conv_*cos(alpha_)*ux               // ux
          + conv_*sin(alpha_)*uy;              // uy
      }
    }

  } else if( RhsType_ == "exact_rhs_recirc_2d" ) {

    // need to set a_ and b_ too
    if( conv_ == UNDEF ) conv_ = 1;
    if( diff_ == UNDEF ) diff_ = 1e-5;

    SetupCartesianGrid2D();

    double hx = lx_/(nx_+1);
    double hy = ly_/(ny_+1);

    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int ix, iy;
      ix = (MyGlobalElements[i])%nx_;
      iy = (MyGlobalElements[i] - ix)/nx_;
      double x = hx*(ix+1);
      double y = hy*(iy+1);
      double u, ux, uy, uxx, uyy;
      ExactSolQuadXY(x,y,u,ux,uy,uxx,uyy);

      for (int j = 0 ; j < NumVectors_ ; ++j) {
        (*rhs_)[j][i] =  -diff_*( uxx + uyy )        // -b_ \nabla u
          + conv_*4*x*(x-1.)*(1.-2*y)*ux          // ux
          - conv_*4*y*(y-1.)*(1.-2*x)*uy;         // uy
      }
    }

  } else if( RhsType_ == "exact_rhs_laplace_2d" ) {

    SetupCartesianGrid2D();

    double hx = lx_/(nx_+1);
    double hy = ly_/(ny_+1);

    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int ix, iy;
      ix = (MyGlobalElements[i])%nx_;
      iy = (MyGlobalElements[i] - ix)/nx_;
      double x = hx*(ix+1);
      double y = hy*(iy+1);
      double u, ux, uy, uxx, uyy;
      ExactSolQuadXY(x,y,u,ux,uy,uxx,uyy);

      for (int j = 0 ; j < NumVectors_ ; ++j) {
        for (int jj = 0 ; jj < NumVectors_ ; ++jj) {
          (*rhs_)[jj][i] = uxx+uyy;
        }
      }
    }

  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        ErrorMsg << "RHS type not correct (" << RhsType_ << ")"
        );
  }

  if( verbose_ ) {
    cout << OutputMsg << "Time to create RHS (matvec): "
      << Time.ElapsedTime() << " (s)\n";
  }

  return;
}

void Trilinos_Util::CrsMatrixGallery::CreateRHS(void)
{
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_)
    TCreateRHS<long long>();
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_)
      TCreateRHS<int>();
    else
#endif
      throw "Trilinos_Util::CrsMatrixGallery::CreateRHS: failed";
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateEye(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `eye'...\n";
  }

  a_ = 1.0;
  CreateMatrixDiag<int_type>();
  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixDiag(void)
{

  // default value if not otherwise specified
  if( a_ == UNDEF ) a_ = 1;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `diag'...\n";
    cout << OutputMsg << "Diagonal element = " << a_ << endl;
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,1);
  double Value;
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_; ++i ) {

    int_type Indices = MyGlobalElements[i];
    Value = a_;

    matrix_->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);

  }

  matrix_->FillComplete();

  return;

}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixTriDiag(void)
{
  // default value if not otherwise specified
  if( a_ == UNDEF ) a_ = 2;
  if( b_ == UNDEF ) b_ = 1;
  if( c_ == UNDEF ) c_ = 1;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `tridiag'...\n";
    cout << OutputMsg << "Row is [" << b_ << ", " << a_ << ", " << c_ << "]\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,3);

  double * Values = new double[2];
  int_type * Indices = new int_type[2];
  int NumEntries;

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_; ++i ) {
    if (MyGlobalElements[i]==0) {
      // off-diagonal for first row
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = c_;
    } else if (MyGlobalElements[i] == NumGlobalElements_-1) {
      // off-diagonal for last row
      Indices[0] = NumGlobalElements_-2;
      NumEntries = 1;
      Values[0] = b_;
    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i]-1;
      Values[0] = b_;
      Indices[1] = MyGlobalElements[i]+1;
      Values[1] = c_;
      NumEntries = 2;
    }
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    Values[0] = a_;
    matrix_->InsertGlobalValues(MyGlobalElements[i], 1, Values, MyGlobalElements+i);
  }

  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  matrix_->FillComplete();

  delete [] Values;
  delete [] Indices;

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace1d(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `laplace_1d'...\n";
  }

  a_ = 2.0;
  b_ = -1.0;
  c_ = -1.0;

  CreateMatrixTriDiag<int_type>();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace1dNeumann(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `laplace_1d_n'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,3);

  double *Values = new double[2];
  int_type *Indices = new int_type[2];
  int NumEntries;

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_; ++i ) {
    if (MyGlobalElements[i]==0) {
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = -1.0;
    } else if (MyGlobalElements[i] == NumGlobalElements_-1) {
      Indices[0] = NumGlobalElements_-2;
      NumEntries = 1;
      Values[0] = -1.0;
    } else {
      Indices[0] = MyGlobalElements[i]-1;
      Values[1] = -1.0;
      Indices[1] = MyGlobalElements[i]+1;
      Values[0] = -1.0;
      NumEntries = 2;
    }
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    if (MyGlobalElements[i]==0 || (MyGlobalElements[i] == NumGlobalElements_-1) )
      Values[0] = 1.0;
    else
      Values[0] = 2.0;

    matrix_->InsertGlobalValues(MyGlobalElements[i], 1, Values, MyGlobalElements+i);
  }

  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  matrix_->FillComplete();

  delete [] Values;
  delete [] Indices;

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil2d(void)
{
  // default values if not otherwise specified

  if( a_ == UNDEF ) a_ = 4;
  if( b_ == UNDEF ) b_ = 1;
  if( c_ == UNDEF ) c_ = 1;
  if( d_ == UNDEF ) d_ = 1;
  if( e_ == UNDEF ) e_ = 1;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
    cout << OutputMsg << "with values: a=" << a_ << ", b=" << b_ << ", c=" << c_
      << ", d=" << d_ << ", e=" << e_ << endl;
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);

  // Add  rows one-at-a-time

  double Values[4], diag;
  int_type Indices[4];

  //    e
  //  b a c
  //    d
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements[i], nx_, ny_,
        left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = b_;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = c_;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d_;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e_;
      ++NumEntries;
    }
    // put the off-diagonal entries
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries,
        Values, Indices);
    // Put in the diagonal entry
    diag = a_;

    matrix_->InsertGlobalValues(MyGlobalElements[i], 1,
        &diag, MyGlobalElements+i);
  }
  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil2dVector(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);

  // Add  rows one-at-a-time

  double Values[4], diag;
  int_type Indices[4];

  //    e
  //  b a c
  //    d
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
  for( int i=0 ; i<NumMyElements_; ++i ) {

    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements[i], nx_, ny_,
        left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries,
        Values, Indices);
    // Put in the diagonal entry
    diag = (*VectorA_)[i];

    matrix_->InsertGlobalValues(MyGlobalElements[i], 1,
        &diag, MyGlobalElements+i);
  }
  matrix_->FillComplete();

  return;

}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace2d_BC(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `laplace_2d_bc'...\n";
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);

  // Add  rows one-at-a-time

  double Values[4], diag;
  int_type Indices[4];

  //    e
  //  b a c
  //    d
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
  for( int i=0 ; i<NumMyElements_; ++i ) {

    //bool isBorder = false;

    //int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements[i], nx_, ny_,
        left, right, lower, upper);

    // any border node gets only diagonal entry
    if ((left == -1) || (right == -1) || (lower == -1) || (upper == -1)) {
      diag = 1.;
    }
    else {

      Indices[0] = left;
      Values[0] = -1.0;

      Indices[1] = right;
      Values[1] = -1.0;

      Indices[2] = lower;
      Values[2] = -1.0;

      Indices[3] = upper;
      Values[3] = -1.0;

      // put the off-diagonal entries
      matrix_->InsertGlobalValues(MyGlobalElements[i], 4, Values, Indices);

      diag = 4.0;
    }

    // new diagonal guy
    matrix_->InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
  }
  matrix_->FillComplete();

  return;

}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace2dNeumann(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `laplace_2d_n'...\n";
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);

  // Add  rows one-at-a-time

  double Values[4], diag;
  int_type Indices[4];

  //    e
  //  b a c
  //    d
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
  for( int i=0 ; i<NumMyElements_; ++i ) {

    bool isBorder = false;

    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements[i], nx_, ny_,
        left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;

    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;

    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;

    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    } else isBorder = true;

    // put the off-diagonal entries
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    if( isBorder ) diag = 2;
    else           diag = 4;

    matrix_->InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
  }
  matrix_->FillComplete();

  return;

}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace2d_9pt(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `laplace_2d_9pt'...\n";
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,9);

  // Add  rows one-at-a-time

  double Values[8], diag;
  for( int i=0 ; i<8 ; ++i ) Values[i] = -1.0;
  int_type Indices[8];

  diag = 8.0;

  //  z3  e  z4
  //   b  a  c
  //  z1  d  z2
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements[i], nx_, ny_,
        left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      ++NumEntries;
    }
    if( left != -1 && lower != -1 ) {
      Indices[NumEntries] = lower-1;
      ++NumEntries;
    }
    if( right != -1 && lower != -1 ) {
      Indices[NumEntries] = lower+1;
      ++NumEntries;
    }
    if( left != -1 && upper != -1 ) {
      Indices[NumEntries] = upper-1;
      ++NumEntries;
    }
    if( right != -1 && upper != -1 ) {
      Indices[NumEntries] = upper+1;
      ++NumEntries;
    }

    // put the off-diagonal entries
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries,
        Values, Indices);
    // Put in the diagonal entry
    matrix_->InsertGlobalValues(MyGlobalElements[i], 1,
        &diag, MyGlobalElements+i);
  }
  matrix_->FillComplete();

  return;

}
// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace2d(void)
{

  SetupCartesianGrid2D();

  double hx = lx_/(nx_+1);
  double hy = ly_/(ny_+1);

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `laplace_2d'...\n";
  }

  if( Scaling ) {
    if( Scaling ) {
      cout << OutputMsg << "hx = " << hx << ", hy = " << hy << endl;
    }
    a_ = 2.0/(hx*hx) + 2.0/(hy*hy);
    b_ = -1.0/(hx*hx);
    c_ = -1.0/(hx*hx);
    d_ = -1.0/(hy*hy);
    e_ = -1.0/(hy*hy);
  } else {
    a_ = 4;
    b_ = -1;
    c_ = -1;
    d_ = -1;
    e_ = -1;
  }

  CreateMatrixCrossStencil2d<int_type>();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixRecirc2d(void)
{

  // default values if not specified otherwise

  if( conv_ == UNDEF ) conv_ = 1;
  if( diff_ == UNDEF ) diff_ = 1e-5;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `recirc_2d'...\n";
    cout << OutputMsg << "with convection = " << conv_ << " and diffusion = " << diff_ << endl;
  }

  SetupCartesianGrid2D();

  if( VectorA_ ) delete VectorA_;
  if( VectorB_ ) delete VectorB_;
  if( VectorC_ ) delete VectorC_;
  if( VectorD_ ) delete VectorD_;
  if( VectorE_ ) delete VectorE_;

  if( VectorA_ == NULL )  VectorA_ = new Epetra_Vector(*map_);
  if( VectorB_ == NULL )  VectorB_ = new Epetra_Vector(*map_);
  if( VectorC_ == NULL )  VectorC_ = new Epetra_Vector(*map_);
  if( VectorD_ == NULL )  VectorD_ = new Epetra_Vector(*map_);
  if( VectorE_ == NULL )  VectorE_ = new Epetra_Vector(*map_);

  assert( VectorA_ != NULL ) ;
  assert( VectorB_ != NULL ) ;
  assert( VectorC_ != NULL ) ;
  assert( VectorD_ != NULL ) ;
  assert( VectorE_ != NULL ) ;

  VectorA_->PutScalar(0.0);
  VectorB_->PutScalar(0.0);
  VectorC_->PutScalar(0.0);
  VectorD_->PutScalar(0.0);
  VectorE_->PutScalar(0.0);

  double hx = lx_/(nx_+1);
  double hy = ly_/(ny_+1);

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int ix, iy;
    ix = (MyGlobalElements[i])%nx_;
    iy = (MyGlobalElements[i] - ix)/nx_;
    double x = hx*(ix+1);
    double y = hy*(iy+1);
    double ConvX = conv_*4*x*(x-1.)*(1.-2*y)/hx;
    double ConvY = -conv_*4*y*(y-1.)*(1.-2*x)/hy;

    // convection part

    if( ConvX<0 ) {
      (*VectorC_)[i] += ConvX;
      (*VectorA_)[i] -= ConvX;
    } else {
      (*VectorB_)[i] -= ConvX;
      (*VectorA_)[i] += ConvX;
    }

    if( ConvY<0 ) {
      (*VectorE_)[i] += ConvY;
      (*VectorA_)[i] -= ConvY;
    } else {
      (*VectorD_)[i] -= ConvY;
      (*VectorA_)[i] += ConvY;
    }

    // add diffusion part
    (*VectorA_)[i] += diff_*2./(hx*hx) + diff_*2./(hy*hy);
    (*VectorB_)[i] -= diff_/(hx*hx);
    (*VectorC_)[i] -= diff_/(hx*hx);
    (*VectorD_)[i] -= diff_/(hy*hy);
    (*VectorE_)[i] -= diff_/(hy*hy);


  }

  CreateMatrixCrossStencil2dVector<int_type>();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixRecirc2dDivFree(void)
{

  // default values if not specified otherwise

  if( conv_ == UNDEF ) conv_ = 1;
  if( diff_ == UNDEF ) diff_ = 1e-5;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `recirc_2d_divfree'...\n";
    cout << OutputMsg << "with convection = " << conv_ << " and diffusion = " << diff_ << endl;
  }


  SetupCartesianGrid2D();

  if( VectorA_ ) delete VectorA_;
  if( VectorB_ ) delete VectorB_;
  if( VectorC_ ) delete VectorC_;
  if( VectorD_ ) delete VectorD_;
  if( VectorE_ ) delete VectorE_;

  if( VectorA_ == NULL )  VectorA_ = new Epetra_Vector(*map_);
  if( VectorB_ == NULL )  VectorB_ = new Epetra_Vector(*map_);
  if( VectorC_ == NULL )  VectorC_ = new Epetra_Vector(*map_);
  if( VectorD_ == NULL )  VectorD_ = new Epetra_Vector(*map_);
  if( VectorE_ == NULL )  VectorE_ = new Epetra_Vector(*map_);

  VectorA_->PutScalar(0.0);
  VectorB_->PutScalar(0.0);
  VectorC_->PutScalar(0.0);
  VectorD_->PutScalar(0.0);
  VectorE_->PutScalar(0.0);

  double hx = lx_/(nx_+1);
  double hy = ly_/(ny_+1);

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int ix, iy;
    ix = (MyGlobalElements[i])%nx_;
    iy = (MyGlobalElements[i] - ix)/nx_;
    double x = hx*(ix+1);
    double y = hy*(iy+1);
    double ConvX = conv_*2*y*(1.-x*x)/hx;
    double ConvY = -conv_*2*x*(1.-y*y)/hy;

    // convection part

    if( ConvX<0 ) {
      (*VectorC_)[i] += ConvX;
      (*VectorA_)[i] -= ConvX;
    } else {
      (*VectorB_)[i] -= ConvX;
      (*VectorA_)[i] += ConvX;
    }

    if( ConvY<0 ) {
      (*VectorE_)[i] += ConvY;
      (*VectorA_)[i] -= ConvY;
    } else {
      (*VectorD_)[i] -= ConvY;
      (*VectorA_)[i] += ConvY;
    }

    // add diffusion part
    (*VectorA_)[i] += diff_*2./(hx*hx) + diff_*2./(hy*hy);
    (*VectorB_)[i] -= diff_/(hx*hx);
    (*VectorC_)[i] -= diff_/(hx*hx);
    (*VectorD_)[i] -= diff_/(hy*hy);
    (*VectorE_)[i] -= diff_/(hy*hy);


  }

  CreateMatrixCrossStencil2d<int_type>();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixUniFlow2d(void)
{

  // default values if not specified otherwise

  if( conv_ == UNDEF  ) conv_ = 1;
  if( diff_ == UNDEF  ) diff_ = 1e-5;
  if( alpha_ == UNDEF ) alpha_ = 0;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `uni_flow_2d'...\n";
    cout << OutputMsg << "with convection = " << conv_ << ", diffusion = " << diff_ << endl;
    cout << OutputMsg << "and alpha = " << alpha_ << endl;

  }

  SetupCartesianGrid2D();

  if( VectorA_ ) delete VectorA_;
  if( VectorB_ ) delete VectorB_;
  if( VectorC_ ) delete VectorC_;
  if( VectorD_ ) delete VectorD_;
  if( VectorE_ ) delete VectorE_;

  if( VectorA_ == NULL )  VectorA_ = new Epetra_Vector(*map_);
  if( VectorB_ == NULL )  VectorB_ = new Epetra_Vector(*map_);
  if( VectorC_ == NULL )  VectorC_ = new Epetra_Vector(*map_);
  if( VectorD_ == NULL )  VectorD_ = new Epetra_Vector(*map_);
  if( VectorE_ == NULL )  VectorE_ = new Epetra_Vector(*map_);

  assert( VectorA_ != NULL ) ;
  assert( VectorB_ != NULL ) ;
  assert( VectorC_ != NULL ) ;
  assert( VectorD_ != NULL ) ;
  assert( VectorE_ != NULL ) ;

  VectorA_->PutScalar(0.0);
  VectorB_->PutScalar(0.0);
  VectorC_->PutScalar(0.0);
  VectorD_->PutScalar(0.0);
  VectorE_->PutScalar(0.0);

  double hx = lx_/(nx_+1);
  double hy = ly_/(ny_+1);

  //int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    //int ix = (MyGlobalElements[i])%nx_;
    //int iy = (MyGlobalElements[i] - ix)/nx_;

    double ConvX = conv_ * cos(alpha_) / hx;
    double ConvY = conv_ * sin(alpha_) / hy;

    // convection part

    if( ConvX<0 ) {
      (*VectorC_)[i] += ConvX;
      (*VectorA_)[i] -= ConvX;
    } else {
      (*VectorB_)[i] -= ConvX;
      (*VectorA_)[i] += ConvX;
    }

    if( ConvY<0 ) {
      (*VectorE_)[i] += ConvY;
      (*VectorA_)[i] -= ConvY;
    } else {
      (*VectorD_)[i] -= ConvY;
      (*VectorA_)[i] += ConvY;
    }

    // add diffusion part
    (*VectorA_)[i] += diff_*2./(hx*hx) + diff_*2./(hy*hy);
    (*VectorB_)[i] -= diff_/(hx*hx);
    (*VectorC_)[i] -= diff_/(hx*hx);
    (*VectorD_)[i] -= diff_/(hy*hy);
    (*VectorE_)[i] -= diff_/(hy*hy);


  }

  CreateMatrixCrossStencil2dVector<int_type>();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLaplace3d(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `laplace_3d'...\n";
  }

  a_ = 6.0;
  b_ = -1.0;
  c_ = -1.0;
  d_ = -1.0;
  e_ = -1.0;
  f_ = -1.0;
  g_ = -1.0;

  CreateMatrixCrossStencil3d<int_type>();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixStretched2d(void)
{

  if( epsilon_ == UNDEF ) epsilon_ = 1e-5;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `stretched_2d'...\n";
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;
  double diag = 8.0;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,9);

  // Add  rows one-at-a-time

  double Values[8];
  int_type Indices[8];

  //  z1  d  z2
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements[i], nx_, ny_,
        left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = 2.0-epsilon_;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = 2.0-epsilon_;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = -4.0+epsilon_;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = -4.0+epsilon_;
      ++NumEntries;
    }
    if( left != -1 && lower != -1 ) {
      Indices[NumEntries] = lower-1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( right != -1 && lower != -1 ) {
      Indices[NumEntries] = lower+1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( left != -1 && upper != -1 ) {
      Indices[NumEntries] = upper-1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( right != -1 && upper != -1 ) {
      Indices[NumEntries] = upper+1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }

    // put the off-diagonal entries
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries,
        Values, Indices);
    // Put in the diagonal entry
    matrix_->InsertGlobalValues(MyGlobalElements[i], 1,
        &diag, MyGlobalElements+i);
  }
  matrix_->FillComplete();

  return;

}
// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil3d(void)
{

  // default values if not otherwise specified

  if( a_ == UNDEF ) a_ = 7;
  if( b_ == UNDEF ) b_ = 1;
  if( c_ == UNDEF ) c_ = 1;
  if( d_ == UNDEF ) d_ = 1;
  if( e_ == UNDEF ) e_ = 1;
  if( f_ == UNDEF ) f_ = 1;
  if( g_ == UNDEF ) g_ = 1;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_3d'...\n";
    cout << OutputMsg << "with values: a=" << a_ << ", b=" << b_ << ", c=" << c_ << endl
      << OutputMsg << "d=" << d_ << ", e=" << e_ << ", f=" << f_
      << ", g=" << g_ << endl;
  }

  // problem size

  SetupCartesianGrid3D();

  int left, right, lower, upper, below, above;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,7);

  // Add  rows one-at-a-time

  double Values[6], diag;
  int_type Indices[6];

  //    e
  //  b a c
  //    d
  // + f below and g above

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian3d(MyGlobalElements[i], nx_, ny_, nz_,
        left, right, lower, upper, below, above);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = b_;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = c_;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d_;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e_;
      ++NumEntries;
    }
    if( below != -1 ) {
      Indices[NumEntries] = below;
      Values[NumEntries] = f_;
      ++NumEntries;
    }
    if( above != -1 ) {
      Indices[NumEntries] = above;
      Values[NumEntries] = g_;
      ++NumEntries;
    }
    // put the off-diagonal entries
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries,
        Values, Indices);
    // Put in the diagonal entry
    diag = a_;

    matrix_->InsertGlobalValues(MyGlobalElements[i], 1,
        &diag, MyGlobalElements+i);
  }

  matrix_->FillComplete();
  return;

}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCrossStencil3dVector(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_3d'...\n";
  }

  // problem size

  if( nx_ == -1 || ny_ == -1 || nz_ == -1 ) {
    nx_ = (int)pow(1.0*NumGlobalElements_,0.333334);
    ny_ = nx_;
    nz_ = nx_;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        nx_ * ny_ *nz_ != NumGlobalElements_,
        ErrorMsg << "The number of global elements must be a perfect cube\n"
        << ErrorMsg << "(now is " << NumGlobalElements_ << ")."
        );
  }

  int left, right, lower, upper, below, above;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,7);

  // Add  rows one-at-a-time

  double Values[6], diag;
  int_type Indices[6];

  //    e
  //  b a c
  //    d
  // + f below and g above

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian3d(MyGlobalElements[i], nx_, ny_, nz_,
        left, right, lower, upper, below, above);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    if( below != -1 ) {
      Indices[NumEntries] = below;
      Values[NumEntries] = (*VectorF_)[i];
      ++NumEntries;
    }
    if( above != -1 ) {
      Indices[NumEntries] = above;
      Values[NumEntries] = (*VectorG_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries,
        Values, Indices);
    // Put in the diagonal entry
    diag = (*VectorA_)[i];

    matrix_->InsertGlobalValues(MyGlobalElements[i], 1,
        &diag, MyGlobalElements+i);
  }

  matrix_->FillComplete();
  return;
}
// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixLehmer(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `lehmer'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int_type iGlobal = MyGlobalElements[i];
    for( int_type jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      if( iGlobal>=jGlobal) Values[jGlobal] = 1.0*(jGlobal+1)/(iGlobal+1);
      else                  Values[jGlobal] = 1.0*(iGlobal+1)/(jGlobal+1);
    }
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumGlobalElements_, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixMinij(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `minij'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int_type iGlobal = MyGlobalElements[i];
    for( int_type jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      if( iGlobal>=jGlobal ) Values[jGlobal] = 1.0*(jGlobal+1);
      else                   Values[jGlobal] = 1.0*(iGlobal+1);
    }
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumGlobalElements_, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixRis(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `ris'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int_type iGlobal = MyGlobalElements[i];
    for( int_type jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      Values[jGlobal] = 0.5/(NumGlobalElements_ -(iGlobal+1)-(jGlobal+1)+1.5);

    }
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumGlobalElements_, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixHilbert(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `hilbert'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int_type iGlobal = MyGlobalElements[i];
    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      Values[j] = 1.0/((iGlobal+1)+(j+1)-1);
    }
    matrix_->InsertGlobalValues(MyGlobalElements[i], NumGlobalElements_, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixJordblock(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `jordblock'...\n";
  }

  // default values is not specified
  if( a_ == UNDEF ) a_ = 0.1;

  // create matrix

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,2);

  int_type Indices[2];
  double Values[2];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = 0;
    if( MyGlobalElements[i] != NumGlobalElements_-1 ) {
      Indices[NumEntries] = MyGlobalElements[i]+1;
      Values[NumEntries] = 1.0;
      NumEntries++;
    }
    // diagonal contribution
    Indices[NumEntries] = MyGlobalElements[i];
    if( VectorA_ != NULL ) Values[NumEntries] = (*VectorA_)[i];
    else                   Values[NumEntries] = a_;
    NumEntries++;

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixCauchy(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `cauchy'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = NumGlobalElements_;
    int_type iGlobal = MyGlobalElements[i];

    for( int_type jGlobal=0 ; jGlobal<NumGlobalElements_ ; ++jGlobal ) {
      Indices[jGlobal] = jGlobal;
      Values[jGlobal] = 1.0/(iGlobal+1+jGlobal+1);
    }

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixFiedler(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `fiedler'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = NumGlobalElements_;
    int_type iGlobal = MyGlobalElements[i];
    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = abs((double)(iGlobal-j));
    }

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixHanowa(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `hanowa'...\n";
  }

  // default values

  if( a_ == UNDEF ) a_ = -1;

  // problem size
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      NumGlobalElements_ % 2 != 0,
      ErrorMsg << "`hanowa' matrix requires a even number of points"
      );

  int_type half = NumGlobalElements_/2;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,2);

  int_type Indices[2];
  double Values[2];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = 2;
    int_type Global = MyGlobalElements[i];
    Indices[0] = Global;
    if( Global < half ) Indices[1] = Global + half;
    else                Indices[1] = Global - half;
    Values[0] = a_;
    // convert from C style to FORTRAN style
    if( Global < half ) Values[1] = (double) -Global - 1;
    else                Values[1] = (double)  (Global - half) + 1;

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixKMS(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `kms'...\n";
  }

  // default values

  if( a_ == UNDEF ) a_ = 0.5;
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;
    int_type iGlobal = MyGlobalElements[i];
    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      // cast to avoid error: call of overloaded pow(double&, long long int) is ambiguous
      Values[j] = pow(a_, abs((double)(iGlobal-j)));
    }

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixParter(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `parter'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;
    int iGlobal = MyGlobalElements[i];
    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = 1.0/(iGlobal-j+0.5);
    }

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }


  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixPei(void)
{

  // default values if not specified otherwise
  a_ = 1.0;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `pei'...\n";
    cout << OutputMsg << "with value a=" << a_ << endl;
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;

    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = 1.0;
      if( MyGlobalElements[i] == j ) Values[j] += a_;
    }

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixOnes(void)
{

  // default values if not specified otherwise
  if( a_ == UNDEF ) a_ = 1.0;

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `ones'...\n";
    cout << OutputMsg << "with value a=" << a_ << endl;
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;

    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      Values[j] = a_;
    }

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::CreateMatrixVander(void)
{

  if( verbose_ ) {
    cout << OutputMsg << "Creating matrix `vander'...\n";
  }

  // create matrix

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int_type * Indices = new int_type[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  for( int i=0 ; i<NumMyElements_ ; ++i ) {

    int NumEntries = NumGlobalElements_;

    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      Indices[j] = j;
      // cast to avoid error: call of overloaded pow(double&, long long int) is ambiguous
      Values[j] = pow((*VectorA_)[i],(double)(NumGlobalElements_-j-1));
    }

    matrix_->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);

  }

  delete [] Indices;
  delete [] Values;

  matrix_->FillComplete();

  return;
}

// ================================================ ====== ==== ==== == =

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::TReadMatrix(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Reading " << name_ << "  matrix `"
      << FileName_ << "'...\n";
  }

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA;
  Epetra_Vector * readx;
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;

  // Call routine to read in problem from file

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_ && sizeof(int_type) == sizeof(long long)) {
    if( name_ == "hb" )
      Trilinos_Util_ReadHb2Epetra64(FileName_.c_str(), *comm_, readMap,
          readA, readx,
          readb, readxexact);
    else if( name_ == "matrix_market" )
      Trilinos_Util_ReadMatrixMarket2Epetra64(FileName_.c_str(), *comm_,
          readMap, readA, readx,
          readb, readxexact );
    else if( name_ == "triples_sym" )
      Trilinos_Util_ReadTriples2Epetra64(FileName_.c_str(), false, *comm_,
          readMap, readA, readx,
          readb, readxexact );
    else if( name_ == "triples_nonsym" )
      Trilinos_Util_ReadTriples2Epetra64(FileName_.c_str(), true, *comm_,
          readMap, readA, readx,
          readb, readxexact );
    else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg << "problem type not correct (" << name_ << ")"
          );
    }
  }
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_ && sizeof(int_type) == sizeof(int)) {
      if( name_ == "hb" )
        Trilinos_Util_ReadHb2Epetra(FileName_.c_str(), *comm_, readMap,
            readA, readx,
            readb, readxexact);
      else if( name_ == "matrix_market" )
        Trilinos_Util_ReadMatrixMarket2Epetra(FileName_.c_str(), *comm_,
            readMap, readA, readx,
            readb, readxexact );
      else if( name_ == "triples_sym" )
        Trilinos_Util_ReadTriples2Epetra(FileName_.c_str(), false, *comm_,
            readMap, readA, readx,
            readb, readxexact );
      else if( name_ == "triples_nonsym" )
        Trilinos_Util_ReadTriples2Epetra(FileName_.c_str(), true, *comm_,
            readMap, readA, readx,
            readb, readxexact );
      else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            true,
            ErrorMsg << "problem type not correct (" << name_ << ")"
            );
      }
    }
    else
#endif
      throw "Trilinos_Util::CrsMatrixGallery::TReadMatrix: Global indices int vs. long long internal error";

  NumGlobalElements_ = readMap->NumGlobalElements64();

  if( map_ != NULL ) delete map_;

  // create map for matrix. Use the normal function CreateMap
  // if the user has not specified "greedy" as map type.
  // In this latter case, form on proc 0 a map corresponding to
  // the greedy algorithm. This is some kind of graph decomposition
  // stuff, only cheaper (and that does not required any external
  // library)

  if( MapType_ == "greedy" ) {

    int * part = new int[NumGlobalElements_];

    if( comm_->MyPID() == 0 ) {

      int NumProcs = comm_->NumProc();
      int * ElementsPerDomain = new int[NumProcs];
      int * count = new int[NumProcs];

      // define how many nodes have to be put on each proc

      int div = NumGlobalElements_/NumProcs;
      int mod = NumGlobalElements_%NumProcs;

      for( int i=0 ; i<NumProcs ; ++i ) {
        count[i] = 0;
        ElementsPerDomain[i] = div;
        if( i<mod ) ElementsPerDomain[i]++;
      }

      for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) {
        part[i] = -1;
      }

      int MaxNnzPerRow = readA->MaxNumEntries();
      TEUCHOS_ASSERT(MaxNnzPerRow != 0);

      int CrsNumEntries;
      int * CrsIndices;
      double * CrsValues;

      // start from row 0, assigned to domain 0
      int RootNode = 0;
      part[0] = 0;
      int CurrentDomain = 0;

      bool ok = true;

      while( ok ) {

        readA->ExtractMyRowView(RootNode,CrsNumEntries,
            CrsValues,CrsIndices);

        ok = false;

        for( int j=0 ; j<CrsNumEntries ; ++j ) {

          if( count[CurrentDomain] == ElementsPerDomain[CurrentDomain] )
            CurrentDomain++;

          if( part[CrsIndices[j]] == -1 ) {
            part[CrsIndices[j]] = CurrentDomain;
            if(!ok) {
              ok = true;
              RootNode = CrsIndices[j];
            }
            count[CurrentDomain]++;
          }
        }

        // check if some -1 nodes are still available
        if(!ok) {
          for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
            if( part[j] == -1 ) {
              RootNode = j;
              ok = true;
              break;
            }
          }
        }

      }

      delete [] ElementsPerDomain;
      delete [] count;

    }

    // now broadcast on all procs. This might be pretty expensive...
    comm_->Broadcast(part,NumGlobalElements_,0);

    for( int_type j=0 ; j<NumGlobalElements_ ; ++j ) {
      if( part[j] == -1 ) {
        cerr << ErrorMsg << "part[" << j << "] = -1 \n";
      }
    }

    // count the elements assigned to this proc
    int NumMyElements = 0;
    for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) {
      if( part[i] == comm_->MyPID() ) NumMyElements++;
    }

    // get the loc2global list
    int_type * MyGlobalElements = new int_type[NumMyElements];
    int count = 0;
    for( int_type i=0 ; i<NumGlobalElements_ ; ++i ) {
      if( part[i] == comm_->MyPID() ) MyGlobalElements[count++] = i;
    }

    map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,
        0,*comm_);

    delete [] MyGlobalElements;
    delete [] part;

  } else {
    CreateMap();
  }

  // Create Exporter to distribute read-in matrix and vectors
  Epetra_Export exporter(*readMap, *map_);
  matrix_ = new Epetra_CrsMatrix(Copy, *map_, 0);
  StartingSolution_ = new Epetra_MultiVector(*map_, NumVectors_);
  rhs_ = new Epetra_MultiVector(*map_, NumVectors_);
  ExactSolution_ = new Epetra_MultiVector(*map_, NumVectors_);
  StartingSolution_->Export(*readx, exporter, Add);
  rhs_->Export(*readb, exporter, Add);
  ExactSolution_->Export(*readxexact, exporter, Add);
  matrix_->Export(*readA, exporter, Add);

  matrix_->FillComplete();

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;

  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  map_->MyGlobalElementsPtr(MyGlobalElementsPtr<int_type>());

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::GetNeighboursCartesian2d( const int i, const int nx, const int ny,
    int & left, int & right,
    int & lower, int & upper)
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 )
    left = -1;
  else
    left = i-1;
  if( ix == nx-1 )
    right = -1;
  else
    right = i+1;
  if( iy == 0 )
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 )
    upper = -1;
  else
    upper = i+nx;

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::GetNeighboursCartesian3d( const int i, const int nx, const int ny, const int nz,
    int & left, int & right, int & lower, int & upper,
    int & below, int & above )
{

  int ixy, iz;
  ixy = i%(nx*ny);

  iz = (i - ixy)/(nx*ny);

  if( iz == 0 )
    below = -1;
  else
    below = i-nx*ny;
  if( iz == nz-1 )
    above = -1;
  else
    above = i+nx*ny;

  GetNeighboursCartesian2d( ixy, nx, ny, left, right, lower, upper);

  if( left != -1 ) left += iz*(nx*ny);
  if( right != -1 ) right += iz*(nx*ny);
  if( lower != -1 ) lower += iz*(nx*ny);
  if( upper != -1 ) upper += iz*(nx*ny);

  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ZeroOutData()
{
  NumGlobalElements_ = -1;
  nx_ = -1;    ny_ = -1;     nz_ = -1;
  mx_ = -1;    mx_ = -1;     mz_ = -1;

  lx_ = 1.0;
  ly_ = 1.0;
  lz_ = 1.0;

  a_ = UNDEF, b_ = UNDEF, c_ = UNDEF, d_ = UNDEF, e_ = UNDEF, f_ = UNDEF, g_ = UNDEF;
  alpha_ = UNDEF;
  beta_  = UNDEF;
  gamma_ = UNDEF;
  delta_ = UNDEF;
  epsilon_ = UNDEF;

  conv_ = UNDEF;
  diff_ = UNDEF;
  source_ = UNDEF;

  VectorA_ = NULL;
  VectorB_ = NULL;
  VectorC_ = NULL;
  VectorD_ = NULL;
  VectorE_ = NULL;
  VectorF_ = NULL;
  VectorG_ = NULL;

  map_ = NULL;
  matrix_ = NULL;
  ExactSolution_ = NULL;
  StartingSolution_ = NULL;
  rhs_ = NULL;

  MapType_ = "linear";
  ContiguousMap_ = true ;
  ExactSolutionType_ = "constant";
  StartingSolutionType_ = "zero";
  ExpandType_ = "zero_off_diagonal";
  RhsType_ = "from_exact_solution";

  NumPDEEqns_= 1;
  NumVectors_ = 1;

  LinearProblem_ = NULL;

}

void Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors()
{
  PrintMatrixAndVectors(cout);
}

void Trilinos_Util::CrsMatrixGallery::PrintMatrixAndVectors(ostream & os)
{
  if( comm_->MyPID() == 0 ) {
    os << "*** MATRIX ***\n";
  }

  os << *matrix_;

  if( comm_->MyPID() == 0 ) {
    os << "*** RHS ***\n";
  }

  os << *rhs_;

  return;
}

namespace Trilinos_Util {

  ostream & operator << (ostream& os,
      const Trilinos_Util::CrsMatrixGallery & G )
  {

    bool verbose = (G.comm_->MyPID() == 0);

    if( verbose ) {

      os << " * Solving problem " << G.name_ << endl;
      os << " * Number of global elements : " << G.NumGlobalElements_ << endl;
      os << " * Type of Map : " << G.MapType_ << endl;
      os << " * Number of PDEs : " << G.NumPDEEqns_ << endl;

      // CRS stuff
      if( G.matrix_ != NULL ) {
        os << " * the matrix has been created " << endl;
        os << " * Matrix->OperatorDomainMap().NumGlobalElements() = "
          << G.matrix_->OperatorDomainMap().NumGlobalElements64() << endl;
      }
      if( G.ExactSolution_ != NULL )
        os << " * an exact solution (" << G.ExactSolutionType_
          << ") has been created " << endl;
      if( G.rhs_ != NULL )
        os << " * the RHS has been created " << endl;
    }

    //  os << " * On proc " << G.comm_->MyPID() << " there are "
    //     << G.NumMyElements_ << " elements" << endl;

    return os;

  }
} //namespace Trilinos_Util

  template<typename int_type>
void Trilinos_Util::CrsMatrixGallery::TGetCartesianCoordinates(double * & x,
    double * & y,
    double * & z)
{

  if( map_ == NULL ) CreateMap();

  double length = 1.0;
  double delta_x, delta_y, delta_z;

  int_type ix;
  int iy, iz;

  // need coordinates for all elements in ColMap(). This is because often the
  // coordiantes of the so-called external nodes (in the old Aztec notation) are required
  int NumMyElements = matrix_->RowMatrixColMap().NumMyElements();
  int_type * MyGlobalElements;
  matrix_->RowMatrixColMap().MyGlobalElementsPtr(MyGlobalElements);

  if( name_ == "diag" || name_ == "tridiag"  ||
      name_ == "laplace_1d" || name_ == "eye" ) {

    delta_x = length/(nx_-1);

    x = new double[NumMyElements];
    assert( x != 0 );

    for( int i=0 ; i<NumMyElements ; ++i ) {

      ix = MyGlobalElements[i];
      x[i] = delta_x * ix;

    }

  } else  if( name_ == "laplace_2d" || name_ == "cross_stencil_2d"
      || name_ == "laplace_2d_bc"
      || name_ == "laplace_2d_9pt" || name_ == "recirc_2d"
      || name_ == "laplace_2d_n" || name_ == "uni_flow_2d"
      || name_ == "stretched_2d" ) {

    delta_x = lx_/(nx_-1);
    delta_y = ly_/(ny_-1);

    // the user has to deallocate these guys
    x =  new double[NumMyElements];
    y =  new double[NumMyElements];
    assert( x != 0 ); assert( y != 0 );

    for( int i=0 ; i<NumMyElements ; ++i ) {

      ix = MyGlobalElements[i]%nx_;
      iy = (MyGlobalElements[i] - ix)/ny_;

      x[i] = delta_x * ix;
      y[i] = delta_y * iy;

    }

  } else if( name_ == "laplace_3d" || name_ == "cross_stencil_3d" ) {

    delta_x = lx_/(nx_-1);
    delta_y = ly_/(ny_-1);
    delta_z = lz_/(nz_-1);

    x =  new double[NumMyElements];
    y =  new double[NumMyElements];
    z =  new double[NumMyElements];
    assert( x != 0 ); assert( y != 0 ); assert( z != 0 );

    for( int i=0 ; i<NumMyElements ; i++ ) {

      int ixy = MyGlobalElements[i]%(nx_*ny_);
      iz = (MyGlobalElements[i] - ixy)/(nx_*ny_);

      ix = ixy%nx_;
      iy = (ixy - ix)/ny_;

      x[i] = delta_x * ix;
      y[i] = delta_y * iy;
      z[i] = delta_z * iz;

    }

  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        ErrorMsg << "You can build Cartesian coordinates\n"
        << ErrorMsg << "only with one of the following problem_type:\n"
        << ErrorMsg << "<diag> / <tridiag> / <laplace_1d> / <eye>\n"
        << ErrorMsg << "<laplace_2d> / <cross_stencil_2d> / <laplace_2d_9pt> / <recirc_2d>\n"
        << ErrorMsg << "<laplace_2d_n> / <uni_flow_n>\n"
        << ErrorMsg << "<laplace_3d> / <cross_stencil_3d> / <stretched_2d>\n"
        );
  }

  return;

}

void Trilinos_Util::CrsMatrixGallery::
GetCartesianCoordinates(
    double * & x,
    double * & y,
    double * & z
    )
{
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_)
    TGetCartesianCoordinates<long long>(x,y,z);
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_)
      TGetCartesianCoordinates<int>(x,y,z);
    else
#endif
      throw "Trilinos_Util::CrsMatrixGallery::GetCartesianCoordinates: failed";
}


#include <iostream>
#include <fstream>

// ================================================ ====== ==== ==== == =
int Trilinos_Util::CrsMatrixGallery::WriteMatrix( const string & FileName, const bool UseSparse )

{

  // create matrix is not done yet

  if( matrix_ == NULL ) CreateMatrix();

  int NumMyRows = matrix_->NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  long long NumGlobalRows; // global dimension of the problem
  long long GlobalRow;  // row in global ordering
  long long NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = matrix_->NumGlobalRows64();
  NumGlobalNonzeros = matrix_->NumGlobalNonzeros64();

  // print out on cout if no filename is provided

  long long IndexBase = matrix_->IndexBase64(); // MATLAB start from 0
  if( IndexBase == 0 ) IndexBase = 1;

  // write on file the dimension of the matrix

  if( comm_->MyPID() == 0 ) {

    ofstream fout(FileName.c_str());

    if( UseSparse ) {
      fout << "A = spalloc(";
      fout << NumGlobalRows << ',' << NumGlobalRows;
      fout << ',' << NumGlobalNonzeros << ");\n";
    } else {
      fout << "A = zeros(";
      fout << NumGlobalRows << ',' << NumGlobalRows << ");\n";
    }

    fout.close();

  }

  for( int Proc=0 ; Proc<comm_->NumProc() ; ++Proc ) {

    if( comm_->MyPID() == Proc ) {

      ofstream fout(FileName.c_str(),std::ios::app);

      fout << "% On proc " << Proc << ": ";
      fout << NumMyRows << " rows and ";
      fout << matrix_->NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {

        GlobalRow = matrix_->GRID64(MyRow);

        NumNzRow = matrix_->NumMyEntries(MyRow);
        double *Values = new double[NumNzRow];
        int *Indices = new int[NumNzRow];

        matrix_->ExtractMyRowCopy(MyRow, NumNzRow,
            NumEntries, Values, Indices);
        // print out the elements with MATLAB syntax
        for( int j=0 ; j<NumEntries ; ++j ) {
          fout << "A(" << GlobalRow  + IndexBase
            << "," << matrix_->GCID64(Indices[j]) + IndexBase
            << ") = " << Values[j] << ";\n";
        }

        delete [] Values;
        delete [] Indices;

      }

      fout.close();

    }
    comm_->Barrier();
  }

  if( comm_->MyPID() == 0 ) {
    ofstream fout(FileName.c_str(),std::ios::app);
    fout << "%End of Matrix Output\n";
    fout.close();
  }

  return true;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::SetupCartesianGrid2D()
{
  // needs a square number of nodes or
  // nx and ny set
  if( nx_ == -1 || ny_ == -1 ) {
    nx_ = (int)sqrt((double)NumGlobalElements_);
    ny_ = nx_;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        nx_ * ny_ != NumGlobalElements_,
        ErrorMsg << "The number of global elements must be a perfect square\n"
        << ErrorMsg << "otherwise set nx and ny.\n"
        << ErrorMsg << "(now NumGlobalElements = " << NumGlobalElements_ << ")\n"
        );
  }
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::SetupCartesianGrid3D()
{
  // needs a cube number of nodes or
  // nx, ny and nz set
  if( nx_ == -1 || ny_ == -1 || nz_ == -1 ) {
    nx_ = (int)pow((double)NumGlobalElements_,0.333334);
    ny_ = nx_;
    nz_ = nx_;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        nx_ * ny_ * nz_ != NumGlobalElements_,
        ErrorMsg << "The number of global elements must be a perfect cube\n"
        << ErrorMsg << "otherwise set nx, ny, and nz.\n"
        << ErrorMsg << "(now NumGlobalElements = " << NumGlobalElements_ << ")\n"
        );
  }
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ExactSolQuadXY(double x, double y,
    double & u)
{
  u = x*(1.-x)*y*(1.-y);
  return;
}

// ================================================ ====== ==== ==== == =
void Trilinos_Util::CrsMatrixGallery::ExactSolQuadXY(double x, double y,
    double & u,
    double & ux, double & uy,
    double & uxx, double & uyy)
{
  u = x*(1.-x)*y*(1.-y);
  ux = (1-2*x)*y*(1.-y);
  uy = x*(1.-x)*(1.-2*y);
  uxx = -2*(x-x*x);
  uyy = -2*(y-y*y);

  return;
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES // CJ: TODO FIXME for long long

Trilinos_Util::VbrMatrixGallery::~VbrMatrixGallery()
{
  // VBR data
  if( VbrLinearProblem_ != NULL ) delete VbrLinearProblem_;
  if( VbrMatrix_ != NULL ) delete VbrMatrix_;
  if( VbrExactSolution_ != NULL ) delete VbrExactSolution_;
  if( VbrStartingSolution_ != NULL ) delete VbrStartingSolution_;
  if( VbrRhs_ != NULL ) delete VbrRhs_;
  if( BlockMap_ != NULL ) delete BlockMap_;

}
// ================================================ ====== ==== ==== == =
Epetra_MultiVector * Trilinos_Util::VbrMatrixGallery::GetVbrRHS(void)
{
  if( VbrRhs_ == NULL ) CreateVbrRHS();
  return VbrRhs_;
}

// ================================================ ====== ==== ==== == =
Epetra_MultiVector * Trilinos_Util::VbrMatrixGallery::GetVbrExactSolution(void)
{
  if( VbrExactSolution_ == NULL ) CreateVbrExactSolution();
  return VbrExactSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_MultiVector * Trilinos_Util::VbrMatrixGallery::GetVbrStartingSolution(void)
{
  if( VbrStartingSolution_ == NULL ) CreateVbrStartingSolution();
  return VbrStartingSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_VbrMatrix * Trilinos_Util::VbrMatrixGallery::GetVbrMatrix(const int NumPDEEqns)
{
  if( NumPDEEqns != NumPDEEqns_ ) {
    if( BlockMap_ != NULL ) {
      delete BlockMap_;
      BlockMap_ = NULL;
    }
    NumPDEEqns_ = NumPDEEqns;

  }

  return( GetVbrMatrix() );
}

// ================================================ ====== ==== ==== == =
Epetra_VbrMatrix * Trilinos_Util::VbrMatrixGallery::GetVbrMatrix(void)
{
  if( VbrMatrix_ == NULL ) CreateVbrMatrix();
  return VbrMatrix_;
}

// ================================================ ====== ==== ==== == =
Epetra_VbrMatrix & Trilinos_Util::VbrMatrixGallery::GetVbrMatrixRef(void)
{
  if( VbrMatrix_ == NULL ) CreateVbrMatrix();
  return *VbrMatrix_;
}

// ================================================ ====== ==== ==== == =

Epetra_LinearProblem * Trilinos_Util::VbrMatrixGallery::GetVbrLinearProblem(void)
{
  // pointers, not really needed
  Epetra_VbrMatrix * A;
  Epetra_MultiVector * RHS;
  Epetra_MultiVector * StartingSolution;

  A = GetVbrMatrix();
  RHS = GetVbrRHS();
  StartingSolution = GetVbrStartingSolution();

  // create linear problem
  if( VbrLinearProblem_ != NULL ) delete VbrLinearProblem_;
  VbrLinearProblem_ = new Epetra_LinearProblem(A,StartingSolution,RHS);

  return VbrLinearProblem_;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::CreateVbrExactSolution(void)
{
  if( verbose_ )
    cout << OutputMsg << "Creating exact solution (VBR)...\n";

  // check if already have one; in this case delete and rebuild
  if( VbrExactSolution_ != NULL ) delete VbrExactSolution_;
  // need an exact solution for Crs first
  if( ExactSolution_ == NULL ) CreateExactSolution();
  // need a block map first
  if( BlockMap_ == NULL ) CreateBlockMap();
  // now we can expand to the Vbr format
  VbrExactSolution_ = new Epetra_MultiVector(*BlockMap_,NumVectors_);
  for (int k = 0 ; k < NumVectors_ ; ++k)
    for( int j=0 ; j<NumMyElements_ ; j++ )
      for( int i=0 ; i<NumPDEEqns_ ; ++i ) {
        (*VbrExactSolution_)[k][j*NumPDEEqns_+i] = (*ExactSolution_)[k][j];
      }

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::CreateVbrStartingSolution(void)
{
  if( verbose_ )
    cout << OutputMsg << "Creating Starting Solution (VBR)...\n";

  if( VbrStartingSolution_ != NULL ) {
    delete VbrStartingSolution_;
    VbrStartingSolution_ = NULL;
  }

  // need a rhs for crs
  if( StartingSolution_ == NULL ) CreateStartingSolution();
  // need a block map based on map_
  if( BlockMap_ == NULL ) CreateBlockMap();
  // now we can expand to the Vbr format
  VbrStartingSolution_ = new Epetra_MultiVector(*BlockMap_,NumVectors_);
  for (int k = 0 ; k < NumVectors_ ; ++k)
    for( int j=0 ; j<NumMyElements_ ; j++ )
      for( int i=0 ; i<NumPDEEqns_ ; ++i ) {
        (*VbrStartingSolution_)[k][j*NumPDEEqns_+i] = (*StartingSolution_)[k][j];
      }

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::CreateVbrRHS(void)
{
  if( verbose_ )
    cout << OutputMsg << "Creating RHS (VBR)...\n";

  if( VbrRhs_ != NULL ) {
    delete VbrRhs_;
    VbrRhs_ = NULL;
  }

  // need a rhs for crs
  if( rhs_ == NULL ) CreateRHS();
  // need a block map based on map_
  if( BlockMap_ == NULL ) CreateBlockMap();
  // require VbrMatrix to be formed first
  if( VbrMatrix_ == NULL ) CreateVbrMatrix();
  // also need an exact solution
  if( VbrExactSolution_ == NULL )  CreateVbrExactSolution();

  VbrRhs_ = new Epetra_MultiVector( *BlockMap_,NumVectors_);
  VbrMatrix_->Multiply(false,*VbrExactSolution_,*VbrRhs_);

  return;
}

// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::VbrMatrixGallery::TCreateVbrMatrix(void)
{
  if( verbose_ )
    cout << OutputMsg << "Creating VBR matrix...\n";

  if( matrix_ == NULL ) CreateMatrix();
  if( BlockMap_ == NULL ) CreateBlockMap();

  int MaxNnzPerRow = matrix_->MaxNumEntries();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      MaxNnzPerRow == 0,
      ErrorMsg << "something went wrong in `CreateMatrix'\n"
      << ErrorMsg << "MaxNnzPerRow == 0"
      );

  // create a VBR matrix based on BlockMap
  VbrMatrix_ = new Epetra_VbrMatrix(Copy, *BlockMap_,MaxNnzPerRow);

  // size of each VBR block
  int MaxBlockSize = MaxBlkSize_*MaxBlkSize_;

  int CrsNumEntries;
  int * CrsIndices;
  double * CrsValues;

  int_type * VbrIndices = new int_type[MaxNnzPerRow];
  double * VbrValues = new double[MaxBlockSize];
  int BlockRows = NumPDEEqns_;

  // cycle over all the local rows.

#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
#endif

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    // extract Crs row
    int ierr = matrix_->ExtractMyRowView(
        i,CrsNumEntries,
        CrsValues,CrsIndices
        );
    TEUCHOS_ASSERT(ierr == 0);

    // matrix_ is in local form. Need global indices
    for( int kk=0 ; kk<CrsNumEntries ; ++kk)
      VbrIndices[kk] = matrix_->GCID64(CrsIndices[kk]);

    // with VBR matrices, we have to insert one block at time.
    // This required two more instructions, one to start this
    // process (BeginInsertGlobalValues), and another one to
    // commit the end of submissions (EndSubmitEntries).

#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
    // get GID of local row
    int_type GlobalNode = MyGlobalElements[i];
    VbrMatrix_->BeginInsertGlobalValues(GlobalNode, CrsNumEntries, VbrIndices);
#else
    // CJ TODO FIXME: Vbr matrices cannot be 64 bit GID based yet.
    throw "Trilinos_Util::VbrMatrixGallery::CreateVbrMatrix: No support for BeginInsertGlobalValues";
#endif

    int ExpandTypeInt;

    if( ExpandType_ == "zero_off_diagonal" ) ExpandTypeInt=0;
    else if( ExpandType_ == "random_off_diagonal" ) ExpandTypeInt=1;
    else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg << "ExpandType not correct (" << ExpandType_ << ")"
          );
    }
    Epetra_Util Util;
    double r = 0.0;

    for( int ii=0 ; ii<CrsNumEntries ; ++ii ) {

      for( int k=0 ; k<BlockRows ; ++k ) { // rows
        for( int h=0 ; h<BlockRows ; ++h ) { // cols
          if( k == h ) VbrValues[k+h*BlockRows] = CrsValues[ii];
          else {
            switch( ExpandTypeInt ) {
              case 0:
                r = 0.0;
                break;
              case 1:
                // get a double between -1 and 1
                r = Util.RandomDouble();
                // scale it so that the sum of the block off-diagonal
                // is not greater than the block diangonal
                r /= (1.5*CrsValues[ii]*BlockRows);
                break;
            }
            VbrValues[k+h*BlockRows] = r;
          }
        }
      }

      VbrMatrix_->SubmitBlockEntry(VbrValues,BlockRows,BlockRows,BlockRows);

    }

    VbrMatrix_->EndSubmitEntries();
  }

  delete [] VbrIndices;
  delete [] VbrValues;

  VbrMatrix_->FillComplete();

  return;
}

void Trilinos_Util::VbrMatrixGallery::CreateVbrMatrix(void)
{
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_)
    TCreateVbrMatrix<long long>();
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_)
      TCreateVbrMatrix<int>();
    else
#endif
      throw "Trilinos_Util::VbrMatrixGallery::CreateVbrMatrix: failed";
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::ComputeResidualVbr(double* residual)
{
  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) CreateVbrRHS();

  Epetra_MultiVector Ax(*BlockMap_,NumVectors_);
  VbrMatrix_->Multiply(false, *VbrStartingSolution_, Ax);
  Ax.Update(1.0, *VbrRhs_, -1.0);
  Ax.Norm2(residual);

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::ComputeDiffBetweenStartingAndExactSolutionsVbr(double* residual)
{
  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) CreateVbrRHS();

  Epetra_MultiVector temp(*BlockMap_,NumVectors_);

  temp.Update(1.0, *VbrExactSolution_, -1.0, *VbrStartingSolution_, 0.0);
  temp.Norm2(residual);

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors(ostream & os)
{
  if( comm_->MyPID() == 0 ) {
    os << "*** MATRIX (VBR) ***\n";
  }

  os << *VbrMatrix_;

  if( comm_->MyPID() == 0 ) {
    os << "*** RHS (VBR) ***\n";
  }

  os << *VbrRhs_;

  return;
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util::VbrMatrixGallery::PrintVbrMatrixAndVectors()
{
  PrintVbrMatrixAndVectors(cout);
}

// ================================================ ====== ==== ==== == =
const Epetra_BlockMap * Trilinos_Util::VbrMatrixGallery::GetBlockMap(void)
{
  if( BlockMap_ == NULL ) CreateBlockMap();

  return BlockMap_;
}

// ================================================ ====== ==== ==== == =
const Epetra_BlockMap & Trilinos_Util::VbrMatrixGallery::GetBlockMapRef(void)
{
  if( BlockMap_ == NULL ) CreateBlockMap();
  return *BlockMap_;
}


// ================================================ ====== ==== ==== == =
  template<typename int_type>
void Trilinos_Util::VbrMatrixGallery::TCreateBlockMap(void)
{
  if( verbose_ ) {
    cout << OutputMsg << "Creating BlockMap...\n";
  }

  if( map_ == NULL ) CreateMap();

  Epetra_Time Time(*comm_);

  if( NumPDEEqns_ <= 0 ) {
    cerr << ErrorMsg << "NumPDEEqns not correct (" << NumPDEEqns_ << "(\n";
    cerr << ErrorMsg << "Set it to 1\n";
    NumPDEEqns_ = 1;
  }

  MaxBlkSize_ = NumPDEEqns_;

  int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();

  BlockMap_ = new Epetra_BlockMap((int_type) NumGlobalElements_,NumMyElements_,
      MyGlobalElements,
      NumPDEEqns_,(int_type) 0,*comm_);

  if( verbose_ ) {
    cout << OutputMsg << "Time to create BlockMap: "
      << Time.ElapsedTime() << " (s)\n";
  }

  return;
}

void Trilinos_Util::VbrMatrixGallery::CreateBlockMap(void)
{
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(UseLongLong_)
    TCreateBlockMap<long long>();
  else
#endif
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(!UseLongLong_)
      TCreateBlockMap<int>();
    else
#endif
      throw "Trilinos_Util::VbrMatrixGallery::CreateBlockMap: failed";
}

#endif
