/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef MLMEX_H
#define MLMEX_H

#include "ml_config.h"
#include "ml_common.h"

#ifdef HAVE_ML_MLAPI
/* MLAPI Headers */
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelAdaptiveSA.h"
#include "MLAPI_DistributedMatrix.h"
#include "MLAPI_Krylov.h"
#endif

/* ML_Epetra Headers */
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"

#include "ml_RefMaxwell.h"

#ifdef HAVE_ML_MATLAB
#include "mex.h"

/* The big classes of data */
class ml_data_pack;
class ml_data_pack{
public:
  ml_data_pack();
  virtual ~ml_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  virtual int setup(int N,int* rowind,int* colptr, double* vals)=0;

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE
  */
  virtual int status()=0;

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     A       - The matrix to solve with (may not be the one the preconditioned was used for)
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  virtual int solve(Teuchos::ParameterList* TPL, Epetra_CrsMatrix * A, double*b, double*x, int &iters)=0;

  /* Gets the stored Epetra_CrsMatrix if we're in Epetra mode */
  virtual Epetra_CrsMatrix * GetMatrix()=0;

  /* Returns the number of rows */
  virtual int NumMyRows()=0;

  /* Returns the number of cols */
  virtual int NumMyCols()=0;
public:
  int id;
  Teuchos::ParameterList *List;
  double operator_complexity;
  ml_data_pack *next;
};

class mlapi_data_pack:public ml_data_pack{
public:
  mlapi_data_pack();
  ~mlapi_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup(int N,int* rowind,int* colptr, double* vals);

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE
  */
  int status();


  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O] (NOT IMPLEMENTED)
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x, int &iters);

  /* Returns the number of rows */
  int NumMyRows(){return A->NumMyRows();}

  /* Returns the number of cols */
  int NumMyCols(){return A->NumMyCols();}

private:
  int solve(Teuchos::ParameterList* TPL, Epetra_CrsMatrix * A, double*b, double*x, int &iters) {return 0;}
  Epetra_CrsMatrix * GetMatrix() {return 0;}

  MLAPI::Space * FineSpace;
  MLAPI::DistributedMatrix * A;
  MLAPI::MultiLevelAdaptiveSA *Prec;
};


class ml_epetra_data_pack:public ml_data_pack{
public:
  ml_epetra_data_pack();
  ~ml_epetra_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup(int N,int* rowind,int* colptr, double* vals);

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE
  */
  int status();

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     A       - The matrix to solve with (may not be the one the preconditioned was used for)
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, Epetra_CrsMatrix * Amat, double*b, double*x,int &iters);

  /* GetPreconditioner - returns a pointer to the preconditioner */
  ML_Epetra::MultiLevelPreconditioner* GetPreconditioner(){return Prec;}

  /*GetMatrix - Returns a pointer to the matrix */
  Epetra_CrsMatrix * GetMatrix(){return A;}

  /* Returns the number of rows */
  int NumMyRows(){return A->NumMyRows();}

  /* Returns the number of cols */
  int NumMyCols(){return A->NumMyCols();}
public:
private:
  Epetra_CrsMatrix * A;
  ML_Epetra::MultiLevelPreconditioner *Prec;
};


class ml_maxwell_data_pack:public ml_data_pack{
public:
  ml_maxwell_data_pack();
  ~ml_maxwell_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     name    - Name of matrix to create [I]
     mxa     - mxArray containing data
     rewrap_ints - if the ints don't line up correctly
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup_matrix(const char * name, const mxArray * mxa, bool rewrap_ints);
  /*
    Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup_preconditioner();

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE
  */
  int status();

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     A       - The matrix to solve with (may not be the one the preconditioned was used for)
     b       - RHS vector [I]
     x       - solution vector [O]
     iters   - number of iterations taken [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, Epetra_CrsMatrix *Amat, double*b, double*x,int &iters);

  /* GetPreconditioner - returns a pointer to the preconditioner */
  ML_Epetra::RefMaxwellPreconditioner* GetPreconditioner(){return Prec;}

  /*GetMatrix - Returns a pointer to the matrix */
  Epetra_CrsMatrix * GetMatrix(){return EdgeMatrix;}

  /* Returns the number of rows */
  int NumMyRows(){return EdgeMatrix->NumMyRows();}

  /* Returns the number of cols */
  int NumMyCols(){return EdgeMatrix->NumMyCols();}

private:
  // Disabled
  int setup(int N,int* rowind,int* colptr, double* vals){return 0;}


  Epetra_CrsMatrix *EdgeMatrix, *GradMatrix, *NodeMatrix, *DummyMatrix, *DummyMatrix2;
  //  ML_Epetra::MultiLevelPreconditioner *Prec;
  ML_Epetra::RefMaxwellPreconditioner *Prec;

};

#endif
#endif
