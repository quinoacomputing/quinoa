
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


#include "CTrilinos_config.h"

#include "CTrilinos_enums.h"
#include "CEpetra_LAPACK.h"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "Epetra_LAPACK.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_LAPACK */
Table<Epetra_LAPACK>& tableOfLAPACKs()
{
    static Table<Epetra_LAPACK> loc_tableOfLAPACKs(CT_Epetra_LAPACK_ID);
    return loc_tableOfLAPACKs;
}


} // namespace


//
// Definitions from CEpetra_LAPACK.h
//


extern "C" {


CT_Epetra_LAPACK_ID_t Epetra_LAPACK_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_LAPACK_Generalize ( 
  CT_Epetra_LAPACK_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id);
}

CT_Epetra_LAPACK_ID_t Epetra_LAPACK_Create (  )
{
    return CEpetra::storeNewLAPACK(new Epetra_LAPACK());
}

CT_Epetra_LAPACK_ID_t Epetra_LAPACK_Duplicate ( 
  CT_Epetra_LAPACK_ID_t LAPACKID )
{
    const Teuchos::RCP<const Epetra_LAPACK> LAPACK = CEpetra::getConstLAPACK(
        LAPACKID);
    return CEpetra::storeNewLAPACK(new Epetra_LAPACK(*LAPACK));
}

void Epetra_LAPACK_Destroy ( CT_Epetra_LAPACK_ID_t * selfID )
{
    CEpetra::removeLAPACK(selfID);
}

void Epetra_LAPACK_POTRF_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  float * A, const int LDA, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POTRF(UPLO, N, A, LDA, INFO);
}

void Epetra_LAPACK_POTRF_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  double * A, const int LDA, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POTRF(UPLO, N, A, LDA, INFO);
}

void Epetra_LAPACK_POTRS_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const int NRHS, const float * A, const int LDA, float * X, 
  const int LDX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POTRS(UPLO, N, NRHS, A, LDA, X, LDX, 
        INFO);
}

void Epetra_LAPACK_POTRS_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const int NRHS, const double * A, const int LDA, double * X, 
  const int LDX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POTRS(UPLO, N, NRHS, A, LDA, X, LDX, 
        INFO);
}

void Epetra_LAPACK_POTRI_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  float * A, const int LDA, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POTRI(UPLO, N, A, LDA, INFO);
}

void Epetra_LAPACK_POTRI_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  double * A, const int LDA, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POTRI(UPLO, N, A, LDA, INFO);
}

void Epetra_LAPACK_POCON_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const float * A, const int LDA, const float ANORM, float * RCOND, 
  float * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POCON(UPLO, N, A, LDA, ANORM, RCOND, WORK, 
        IWORK, INFO);
}

void Epetra_LAPACK_POCON_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const double * A, const int LDA, const double ANORM, 
  double * RCOND, double * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POCON(UPLO, N, A, LDA, ANORM, RCOND, WORK, 
        IWORK, INFO);
}

void Epetra_LAPACK_POSV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const int NRHS, float * A, const int LDA, float * X, 
  const int LDX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POSV(UPLO, N, NRHS, A, LDA, X, LDX, INFO);
}

void Epetra_LAPACK_POSV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const int NRHS, double * A, const int LDA, double * X, 
  const int LDX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POSV(UPLO, N, NRHS, A, LDA, X, LDX, INFO);
}

void Epetra_LAPACK_POEQU_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const float * A, 
  const int LDA, float * S, float * SCOND, float * AMAX, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POEQU(N, A, LDA, S, SCOND, AMAX, INFO);
}

void Epetra_LAPACK_POEQU_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const double * A, 
  const int LDA, double * S, double * SCOND, double * AMAX, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POEQU(N, A, LDA, S, SCOND, AMAX, INFO);
}

void Epetra_LAPACK_PORFS_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const int NRHS, const float * A, const int LDA, const float * AF, 
  const int LDAF, const float * B, const int LDB, float * X, 
  const int LDX, float * FERR, float * BERR, float * WORK, 
  int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->PORFS(UPLO, N, NRHS, A, LDA, AF, LDAF, B, 
        LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO);
}

void Epetra_LAPACK_PORFS_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  const int NRHS, const double * A, const int LDA, 
  const double * AF, const int LDAF, const double * B, 
  const int LDB, double * X, const int LDX, double * FERR, 
  double * BERR, double * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->PORFS(UPLO, N, NRHS, A, LDA, AF, LDAF, B, 
        LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO);
}

void Epetra_LAPACK_POSVX_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char FACT, const char UPLO, 
  const int N, const int NRHS, float * A, const int LDA, float * AF, 
  const int LDAF, const char EQUED, float * S, float * B, 
  const int LDB, float * X, const int LDX, float * RCOND, 
  float * FERR, float * BERR, float * WORK, int * IWORK, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POSVX(FACT, UPLO, N, NRHS, A, LDA, AF, 
        LDAF, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
}

void Epetra_LAPACK_POSVX_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char FACT, const char UPLO, 
  const int N, const int NRHS, double * A, const int LDA, 
  double * AF, const int LDAF, const char EQUED, double * S, 
  double * B, const int LDB, double * X, const int LDX, 
  double * RCOND, double * FERR, double * BERR, double * WORK, 
  int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->POSVX(FACT, UPLO, N, NRHS, A, LDA, AF, 
        LDAF, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO);
}

void Epetra_LAPACK_GELS_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int M, 
  const int N, const int NRHS, double * A, const int LDA, 
  double * B, const int LDB, double * WORK, const int LWORK, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GELS(TRANS, M, N, NRHS, A, LDA, B, LDB, 
        WORK, LWORK, INFO);
}

void Epetra_LAPACK_GETRF_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, float * A, 
  const int LDA, int * IPIV, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GETRF(M, N, A, LDA, IPIV, INFO);
}

void Epetra_LAPACK_GETRF_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  double * A, const int LDA, int * IPIV, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GETRF(M, N, A, LDA, IPIV, INFO);
}

void Epetra_LAPACK_GEQRF_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, float * A, 
  const int LDA, float * TAU, float * WORK, const int lwork, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEQRF(M, N, A, LDA, TAU, WORK, lwork, 
        INFO);
}

void Epetra_LAPACK_GEQRF_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  double * A, const int LDA, double * TAU, double * WORK, 
  const int lwork, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEQRF(M, N, A, LDA, TAU, WORK, lwork, 
        INFO);
}

void Epetra_LAPACK_GETRS_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  const int NRHS, const float * A, const int LDA, const int * IPIV, 
  float * X, const int LDX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GETRS(TRANS, N, NRHS, A, LDA, IPIV, X, 
        LDX, INFO);
}

void Epetra_LAPACK_GETRS_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  const int NRHS, const double * A, const int LDA, const int * IPIV, 
  double * X, const int LDX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GETRS(TRANS, N, NRHS, A, LDA, IPIV, X, 
        LDX, INFO);
}

void Epetra_LAPACK_GETRI_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, float * A, 
  const int LDA, int * IPIV, float * WORK, const int * LWORK, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GETRI(N, A, LDA, IPIV, WORK, LWORK, INFO);
}

void Epetra_LAPACK_GETRI_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, double * A, 
  const int LDA, int * IPIV, double * WORK, const int * LWORK, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GETRI(N, A, LDA, IPIV, WORK, LWORK, INFO);
}

void Epetra_LAPACK_GECON_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char NORM, const int N, 
  const float * A, const int LDA, const float ANORM, float * RCOND, 
  float * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GECON(NORM, N, A, LDA, ANORM, RCOND, WORK, 
        IWORK, INFO);
}

void Epetra_LAPACK_GECON_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char NORM, const int N, 
  const double * A, const int LDA, const double ANORM, 
  double * RCOND, double * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GECON(NORM, N, A, LDA, ANORM, RCOND, WORK, 
        IWORK, INFO);
}

void Epetra_LAPACK_GESV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const int NRHS, 
  float * A, const int LDA, int * IPIV, float * X, const int LDX, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESV(N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}

void Epetra_LAPACK_GESV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const int NRHS, 
  double * A, const int LDA, int * IPIV, double * X, const int LDX, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESV(N, NRHS, A, LDA, IPIV, X, LDX, INFO);
}

void Epetra_LAPACK_GEEQU_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  const float * A, const int LDA, float * R, float * C, 
  float * ROWCND, float * COLCND, float * AMAX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEEQU(M, N, A, LDA, R, C, ROWCND, COLCND, 
        AMAX, INFO);
}

void Epetra_LAPACK_GEEQU_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  const double * A, const int LDA, double * R, double * C, 
  double * ROWCND, double * COLCND, double * AMAX, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEEQU(M, N, A, LDA, R, C, ROWCND, COLCND, 
        AMAX, INFO);
}

void Epetra_LAPACK_GERFS_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  const int NRHS, const float * A, const int LDA, const float * AF, 
  const int LDAF, const int * IPIV, const float * B, const int LDB, 
  float * X, const int LDX, float * FERR, float * BERR, 
  float * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GERFS(TRANS, N, NRHS, A, LDA, AF, LDAF, 
        IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO);
}

void Epetra_LAPACK_GERFS_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  const int NRHS, const double * A, const int LDA, 
  const double * AF, const int LDAF, const int * IPIV, 
  const double * B, const int LDB, double * X, const int LDX, 
  double * FERR, double * BERR, double * WORK, int * IWORK, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GERFS(TRANS, N, NRHS, A, LDA, AF, LDAF, 
        IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO);
}

void Epetra_LAPACK_GESVX_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char FACT, const char TRANS, 
  const int N, const int NRHS, float * A, const int LDA, float * AF, 
  const int LDAF, int * IPIV, const char EQUED, float * R, 
  float * C, float * B, const int LDB, float * X, const int LDX, 
  float * RCOND, float * FERR, float * BERR, float * WORK, 
  int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESVX(FACT, TRANS, N, NRHS, A, LDA, AF, 
        LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, 
        IWORK, INFO);
}

void Epetra_LAPACK_GESVX_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char FACT, const char TRANS, 
  const int N, const int NRHS, double * A, const int LDA, 
  double * AF, const int LDAF, int * IPIV, const char EQUED, 
  double * R, double * C, double * B, const int LDB, double * X, 
  const int LDX, double * RCOND, double * FERR, double * BERR, 
  double * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESVX(FACT, TRANS, N, NRHS, A, LDA, AF, 
        LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, 
        IWORK, INFO);
}

void Epetra_LAPACK_GEHRD_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  const int IHI, float * A, const int LDA, float * TAU, 
  float * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEHRD(N, ILO, IHI, A, LDA, TAU, WORK, 
        LWORK, INFO);
}

void Epetra_LAPACK_GEHRD_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  const int IHI, double * A, const int LDA, double * TAU, 
  double * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEHRD(N, ILO, IHI, A, LDA, TAU, WORK, 
        LWORK, INFO);
}

void Epetra_LAPACK_HSEQR_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOB, const char COMPZ, 
  const int N, const int ILO, const int IHI, float * H, 
  const int LDH, float * WR, float * WI, float * Z, const int LDZ, 
  float * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->HSEQR(JOB, COMPZ, N, ILO, IHI, H, LDH, WR, 
        WI, Z, LDZ, WORK, LWORK, INFO);
}

void Epetra_LAPACK_HSEQR_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOB, const char COMPZ, 
  const int N, const int ILO, const int IHI, double * H, 
  const int LDH, double * WR, double * WI, double * Z, 
  const int LDZ, double * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->HSEQR(JOB, COMPZ, N, ILO, IHI, H, LDH, WR, 
        WI, Z, LDZ, WORK, LWORK, INFO);
}

void Epetra_LAPACK_ORGQR_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  const int K, float * A, const int LDA, float * TAU, float * WORK, 
  const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->ORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, 
        INFO);
}

void Epetra_LAPACK_ORGQR_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  const int K, double * A, const int LDA, double * TAU, 
  double * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->ORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, 
        INFO);
}

void Epetra_LAPACK_ORGHR_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  const int IHI, float * A, const int LDA, float * TAU, 
  float * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->ORGHR(N, ILO, IHI, A, LDA, TAU, WORK, 
        LWORK, INFO);
}

void Epetra_LAPACK_ORGHR_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  const int IHI, double * A, const int LDA, double * TAU, 
  double * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->ORGHR(N, ILO, IHI, A, LDA, TAU, WORK, 
        LWORK, INFO);
}

void Epetra_LAPACK_ORMHR_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char SIDE, const char TRANS, 
  const int M, const int N, const int ILO, const int IHI, 
  const float * A, const int LDA, const float * TAU, float * C, 
  const int LDC, float * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->ORMHR(SIDE, TRANS, M, N, ILO, IHI, A, LDA, 
        TAU, C, LDC, WORK, LWORK, INFO);
}

void Epetra_LAPACK_ORMHR_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char SIDE, const char TRANS, 
  const int M, const int N, const int ILO, const int IHI, 
  const double * A, const int LDA, const double * TAU, double * C, 
  const int LDC, double * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->ORMHR(SIDE, TRANS, M, N, ILO, IHI, A, LDA, 
        TAU, C, LDC, WORK, LWORK, INFO);
}

void Epetra_LAPACK_LARFT_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char DIRECT, 
  const char STOREV, const int N, const int K, double * V, 
  const int LDV, double * TAU, double * T, const int LDT )
{
    CEpetra::getConstLAPACK(selfID)->LARFT(DIRECT, STOREV, N, K, V, LDV, TAU, 
        T, LDT);
}

void Epetra_LAPACK_LARFT_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char DIRECT, 
  const char STOREV, const int N, const int K, float * V, 
  const int LDV, float * TAU, float * T, const int LDT )
{
    CEpetra::getConstLAPACK(selfID)->LARFT(DIRECT, STOREV, N, K, V, LDV, TAU, 
        T, LDT);
}

void Epetra_LAPACK_TREVC_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char SIDE, const char HOWMNY, 
  int * SELECT, const int N, const float * T, const int LDT, 
  float * VL, const int LDVL, float * VR, const int LDVR, 
  const int MM, int * M, float * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->TREVC(SIDE, HOWMNY, SELECT, N, T, LDT, VL, 
        LDVL, VR, LDVR, MM, M, WORK, INFO);
}

void Epetra_LAPACK_TREVC_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char SIDE, const char HOWMNY, 
  int * SELECT, const int N, const double * T, const int LDT, 
  double * VL, const int LDVL, double * VR, const int LDVR, 
  const int MM, int * M, double * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->TREVC(SIDE, HOWMNY, SELECT, N, T, LDT, VL, 
        LDVL, VR, LDVR, MM, M, WORK, INFO);
}

void Epetra_LAPACK_TREXC_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char COMPQ, const int N, 
  float * T, const int LDT, float * Q, const int LDQ, int IFST, 
  int ILST, float * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->TREXC(COMPQ, N, T, LDT, Q, LDQ, IFST, 
        ILST, WORK, INFO);
}

void Epetra_LAPACK_TREXC_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char COMPQ, const int N, 
  double * T, const int LDT, double * Q, const int LDQ, int IFST, 
  int ILST, double * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->TREXC(COMPQ, N, T, LDT, Q, LDQ, IFST, 
        ILST, WORK, INFO);
}

void Epetra_LAPACK_GESVD_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBU, const char JOBVT, 
  const int M, const int N, float * A, const int LDA, float * S, 
  float * U, const int LDU, float * VT, const int LDVT, 
  float * WORK, const int * LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESVD(JOBU, JOBVT, M, N, A, LDA, S, U, 
        LDU, VT, LDVT, WORK, LWORK, INFO);
}

void Epetra_LAPACK_GESVD_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBU, const char JOBVT, 
  const int M, const int N, double * A, const int LDA, double * S, 
  double * U, const int LDU, double * VT, const int LDVT, 
  double * WORK, const int * LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESVD(JOBU, JOBVT, M, N, A, LDA, S, U, 
        LDU, VT, LDVT, WORK, LWORK, INFO);
}

void Epetra_LAPACK_GGSVD_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBU, const char JOBV, 
  const char JOBQ, const int M, const int N, const int P, int * K, 
  int * L, double * A, const int LDA, double * B, const int LDB, 
  double * ALPHA, double * BETA, double * U, const int LDU, 
  double * V, const int LDV, double * Q, const int LDQ, 
  double * WORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GGSVD(JOBU, JOBV, JOBQ, M, N, P, K, L, A, 
        LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK, INFO);
}

void Epetra_LAPACK_GGSVD_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBU, const char JOBV, 
  const char JOBQ, const int M, const int N, const int P, int * K, 
  int * L, float * A, const int LDA, float * B, const int LDB, 
  float * ALPHA, float * BETA, float * U, const int LDU, float * V, 
  const int LDV, float * Q, const int LDQ, float * WORK, 
  int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GGSVD(JOBU, JOBV, JOBQ, M, N, P, K, L, A, 
        LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK, INFO);
}

void Epetra_LAPACK_GEEV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, const char JOBVR, 
  const int N, double * A, const int LDA, double * WR, double * WI, 
  double * VL, const int LDVL, double * VR, const int LDVR, 
  double * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, 
        LDVL, VR, LDVR, WORK, LWORK, INFO);
}

void Epetra_LAPACK_GEEV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, const char JOBVR, 
  const int N, float * A, const int LDA, float * WR, float * WI, 
  float * VL, const int LDVL, float * VR, const int LDVR, 
  float * WORK, const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, 
        LDVL, VR, LDVR, WORK, LWORK, INFO);
}

void Epetra_LAPACK_SPEV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char UPLO, 
  const int N, double * AP, double * W, double * Z, int LDZ, 
  double * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SPEV(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, 
        INFO);
}

void Epetra_LAPACK_SPEV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char UPLO, 
  const int N, float * AP, float * W, float * Z, int LDZ, 
  float * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SPEV(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, 
        INFO);
}

void Epetra_LAPACK_SPGV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, const char JOBZ, 
  const char UPLO, const int N, double * AP, double * BP, 
  double * W, double * Z, const int LDZ, double * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SPGV(ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, 
        LDZ, WORK, INFO);
}

void Epetra_LAPACK_SPGV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, const char JOBZ, 
  const char UPLO, const int N, float * AP, float * BP, float * W, 
  float * Z, const int LDZ, float * WORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SPGV(ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, 
        LDZ, WORK, INFO);
}

void Epetra_LAPACK_SYEV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char UPLO, 
  const int N, double * A, const int LDA, double * W, double * WORK, 
  const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEV(JOBZ, UPLO, N, A, LDA, W, WORK, 
        LWORK, INFO);
}

void Epetra_LAPACK_SYEV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char UPLO, 
  const int N, float * A, const int LDA, float * W, float * WORK, 
  const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEV(JOBZ, UPLO, N, A, LDA, W, WORK, 
        LWORK, INFO);
}

void Epetra_LAPACK_SYEVD_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char UPLO, 
  const int N, double * A, const int LDA, double * W, double * WORK, 
  const int LWORK, int * IWORK, const int LIWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEVD(JOBZ, UPLO, N, A, LDA, W, WORK, 
        LWORK, IWORK, LIWORK, INFO);
}

void Epetra_LAPACK_SYEVD_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char UPLO, 
  const int N, float * A, const int LDA, float * W, float * WORK, 
  const int LWORK, int * IWORK, const int LIWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEVD(JOBZ, UPLO, N, A, LDA, W, WORK, 
        LWORK, IWORK, LIWORK, INFO);
}

void Epetra_LAPACK_SYEVX_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char RANGE, 
  const char UPLO, const int N, double * A, const int LDA, 
  const double * VL, const double * VU, const int * IL, 
  const int * IU, const double ABSTOL, int * M, double * W, 
  double * Z, const int LDZ, double * WORK, const int LWORK, 
  int * IWORK, int * IFAIL, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEVX(JOBZ, RANGE, UPLO, N, A, LDA, VL, 
        VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO);
}

void Epetra_LAPACK_SYEVX_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char RANGE, 
  const char UPLO, const int N, float * A, const int LDA, 
  const float * VL, const float * VU, const int * IL, 
  const int * IU, const float ABSTOL, int * M, float * W, float * Z, 
  const int LDZ, float * WORK, const int LWORK, int * IWORK, 
  int * IFAIL, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEVX(JOBZ, RANGE, UPLO, N, A, LDA, VL, 
        VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO);
}

void Epetra_LAPACK_SYGV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, const char JOBZ, 
  const char UPLO, const int N, double * A, const int LDA, 
  double * B, const int LDB, double * W, double * WORK, 
  const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYGV(ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, 
        W, WORK, LWORK, INFO);
}

void Epetra_LAPACK_SYGV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, const char JOBZ, 
  const char UPLO, const int N, float * A, const int LDA, float * B, 
  const int LDB, float * W, float * WORK, const int LWORK, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYGV(ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, 
        W, WORK, LWORK, INFO);
}

void Epetra_LAPACK_SYGVX_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, const char JOBZ, 
  const char RANGE, const char UPLO, const int N, double * A, 
  const int LDA, double * B, const int LDB, const double * VL, 
  const double * VU, const int * IL, const int * IU, 
  const double ABSTOL, int * M, double * W, double * Z, 
  const int LDZ, double * WORK, const int LWORK, int * IWORK, 
  int * IFAIL, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYGVX(ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, 
        B, LDB, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, 
        IFAIL, INFO);
}

void Epetra_LAPACK_SYGVX_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, const char JOBZ, 
  const char RANGE, const char UPLO, const int N, float * A, 
  const int LDA, float * B, const int LDB, const float * VL, 
  const float * VU, const int * IL, const int * IU, 
  const float ABSTOL, int * M, float * W, float * Z, const int LDZ, 
  float * WORK, const int LWORK, int * IWORK, int * IFAIL, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYGVX(ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, 
        B, LDB, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, 
        IFAIL, INFO);
}

void Epetra_LAPACK_SYEVR_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char RANGE, 
  const char UPLO, const int N, double * A, const int LDA, 
  const double * VL, const double * VU, const int * IL, 
  const int * IU, const double ABSTOL, int * M, double * W, 
  double * Z, const int LDZ, int * ISUPPZ, double * WORK, 
  const int LWORK, int * IWORK, const int LIWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEVR(JOBZ, RANGE, UPLO, N, A, LDA, VL, 
        VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, 
        INFO);
}

void Epetra_LAPACK_SYEVR_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const char RANGE, 
  const char UPLO, const int N, float * A, const int LDA, 
  const float * VL, const float * VU, const int * IL, 
  const int * IU, const float ABSTOL, int * M, float * W, float * Z, 
  const int LDZ, int * ISUPPZ, float * WORK, const int LWORK, 
  int * IWORK, const int LIWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->SYEVR(JOBZ, RANGE, UPLO, N, A, LDA, VL, 
        VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, 
        INFO);
}

void Epetra_LAPACK_GEEVX_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char BALANC, const char JOBVL, 
  const char JOBVR, const char SENSE, const int N, double * A, 
  const int LDA, double * WR, double * WI, double * VL, 
  const int LDVL, double * VR, const int LDVR, int * ILO, int * IHI, 
  double * SCALE, double * ABNRM, double * RCONDE, double * RCONDV, 
  double * WORK, const int LWORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEEVX(BALANC, JOBVL, JOBVR, SENSE, N, A, 
        LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, 
        RCONDV, WORK, LWORK, IWORK, INFO);
}

void Epetra_LAPACK_GEEVX_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char BALANC, const char JOBVL, 
  const char JOBVR, const char SENSE, const int N, float * A, 
  const int LDA, float * WR, float * WI, float * VL, const int LDVL, 
  float * VR, const int LDVR, int * ILO, int * IHI, float * SCALE, 
  float * ABNRM, float * RCONDE, float * RCONDV, float * WORK, 
  const int LWORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GEEVX(BALANC, JOBVL, JOBVR, SENSE, N, A, 
        LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, 
        RCONDV, WORK, LWORK, IWORK, INFO);
}

void Epetra_LAPACK_GESDD_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const int M, 
  const int N, double * A, const int LDA, double * S, double * U, 
  const int LDU, double * VT, const int LDVT, double * WORK, 
  const int LWORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, 
        LDVT, WORK, LWORK, IWORK, INFO);
}

void Epetra_LAPACK_GESDD_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const int M, 
  const int N, float * A, const int LDA, float * S, float * U, 
  const int LDU, float * VT, const int LDVT, float * WORK, 
  const int LWORK, int * IWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GESDD(JOBZ, M, N, A, LDA, S, U, LDU, VT, 
        LDVT, WORK, LWORK, IWORK, INFO);
}

void Epetra_LAPACK_GGEV_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, const char JOBVR, 
  const int N, double * A, const int LDA, double * B, const int LDB, 
  double * ALPHAR, double * ALPHAI, double * BETA, double * VL, 
  const int LDVL, double * VR, const int LDVR, double * WORK, 
  const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GGEV(JOBVL, JOBVR, N, A, LDA, B, LDB, 
        ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO);
}

void Epetra_LAPACK_GGEV_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, const char JOBVR, 
  const int N, float * A, const int LDA, float * B, const int LDB, 
  float * ALPHAR, float * ALPHAI, float * BETA, float * VL, 
  const int LDVL, float * VR, const int LDVR, float * WORK, 
  const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GGEV(JOBVL, JOBVR, N, A, LDA, B, LDB, 
        ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO);
}

void Epetra_LAPACK_GGLSE_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  const int P, double * A, const int LDA, double * B, const int LDB, 
  double * C, double * D, double * X, double * WORK, 
  const int LWORK, int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GGLSE(M, N, P, A, LDA, B, LDB, C, D, X, 
        WORK, LWORK, INFO);
}

void Epetra_LAPACK_GGLSE_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  const int P, float * A, const int LDA, float * B, const int LDB, 
  float * C, float * D, float * X, float * WORK, const int LWORK, 
  int * INFO )
{
    CEpetra::getConstLAPACK(selfID)->GGLSE(M, N, P, A, LDA, B, LDB, C, D, X, 
        WORK, LWORK, INFO);
}

void Epetra_LAPACK_LAMCH_float ( 
  CT_Epetra_LAPACK_ID_t selfID, const char CMACH, float * T )
{
    CEpetra::getConstLAPACK(selfID)->LAMCH(CMACH, *T);
}

void Epetra_LAPACK_LAMCH_double ( 
  CT_Epetra_LAPACK_ID_t selfID, const char CMACH, double * T )
{
    CEpetra::getConstLAPACK(selfID)->LAMCH(CMACH, *T);
}


} // extern "C"


//
// Definitions from CEpetra_LAPACK_Cpp.hpp
//


/* get Epetra_LAPACK from non-const table using CT_Epetra_LAPACK_ID */
const Teuchos::RCP<Epetra_LAPACK>
CEpetra::getLAPACK( CT_Epetra_LAPACK_ID_t id )
{
    if (tableOfLAPACKs().isType(id.table))
        return tableOfLAPACKs().get<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
}

/* get Epetra_LAPACK from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_LAPACK>
CEpetra::getLAPACK( CTrilinos_Universal_ID_t id )
{
    if (tableOfLAPACKs().isType(id.table))
        return tableOfLAPACKs().get<Epetra_LAPACK>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_LAPACK>(id);
}

/* get const Epetra_LAPACK from either the const or non-const table
 * using CT_Epetra_LAPACK_ID */
const Teuchos::RCP<const Epetra_LAPACK>
CEpetra::getConstLAPACK( CT_Epetra_LAPACK_ID_t id )
{
    if (tableOfLAPACKs().isType(id.table))
        return tableOfLAPACKs().getConst<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
}

/* get const Epetra_LAPACK from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_LAPACK>
CEpetra::getConstLAPACK( CTrilinos_Universal_ID_t id )
{
    if (tableOfLAPACKs().isType(id.table))
        return tableOfLAPACKs().getConst<Epetra_LAPACK>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_LAPACK>(id);
}

/* store Epetra_LAPACK (owned) in non-const table */
CT_Epetra_LAPACK_ID_t
CEpetra::storeNewLAPACK( Epetra_LAPACK *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        tableOfLAPACKs().store<Epetra_LAPACK>(pobj, true));
}

/* store Epetra_LAPACK in non-const table */
CT_Epetra_LAPACK_ID_t
CEpetra::storeLAPACK( Epetra_LAPACK *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        tableOfLAPACKs().store<Epetra_LAPACK>(pobj, false));
}

/* store const Epetra_LAPACK in const table */
CT_Epetra_LAPACK_ID_t
CEpetra::storeConstLAPACK( const Epetra_LAPACK *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        tableOfLAPACKs().store<Epetra_LAPACK>(pobj, false));
}

/* remove Epetra_LAPACK from table using CT_Epetra_LAPACK_ID */
void
CEpetra::removeLAPACK( CT_Epetra_LAPACK_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(*id);
    if (tableOfLAPACKs().isType(aid.table))
        tableOfLAPACKs().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(aid);
}

/* remove Epetra_LAPACK from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeLAPACK( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfLAPACKs().isType(aid->table))
        tableOfLAPACKs().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_LAPACK table */
void
CEpetra::purgeLAPACK(  )
{
    tableOfLAPACKs().purge();
}

/* store Epetra_LAPACK in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasLAPACK( const Teuchos::RCP< Epetra_LAPACK > & robj )
{
    return tableOfLAPACKs().alias(robj);
}

/* store const Epetra_LAPACK in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstLAPACK( const Teuchos::RCP< const Epetra_LAPACK > & robj )
{
    return tableOfLAPACKs().alias(robj);
}



