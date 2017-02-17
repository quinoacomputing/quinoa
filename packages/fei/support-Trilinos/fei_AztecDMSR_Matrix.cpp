/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov)
//
// ************************************************************************
// @HEADER
*/


#include <fei_trilinos_macros.hpp>
#include <fei_iostream.hpp>

#ifdef HAVE_FEI_AZTECOO

#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <fei_mpi.h>

#ifndef FEI_SER

#define AZTEC_MPI AZTEC_MPI
#define AZ_MPI AZ_MPI
#ifndef MPI
#define MPI MPI
#endif

#endif

#include <az_aztec.h>

#include <fei_ArrayUtils.hpp>

#include <fei_Aztec_Map.hpp>
#include <fei_Aztec_LSVector.hpp>
#include <fei_AztecDMSR_Matrix.hpp>

#define ADMSR_LOCAL_ROW_ALLOC_LEN(localRow) 1+bindx[localRow+1]-bindx[localRow]

#define ADMSR_LOCAL_ROW_LEN(localRow) 1+rowLengths_[localRow]

namespace fei_trilinos {

  //==============================================================================
  AztecDMSR_Matrix::AztecDMSR_Matrix(fei::SharedPtr<Aztec_Map> map)
    : isFilled_(false),
    isAllocated_(false),
    localOffset_(map->localOffset()),
    localSize_(map->localSize()),
    amap_(map),
    Amat_(NULL),
    arraysAllocated_(false),
    val(NULL),
    bindx(NULL),
    rowLengths_(NULL),
    nnzeros_(0),
    N_update_(map->localSize()),
    tmp_array_(0),
    tmp_array_len_(0),
    dtmp_array_(0),
    dtmp_array_len_(0),
    azTransformed_(false)
  {
    if (N_update_ > 0) {
      rowLengths_ = new int[N_update_];
      for(int i=0; i<N_update_; i++) {
        rowLengths_[i] = 0;
      }
    }

    Amat_ = AZ_matrix_create(N_update_);
  }

  //==============================================================================
  AztecDMSR_Matrix::AztecDMSR_Matrix(const AztecDMSR_Matrix& src)
    : isFilled_(src.isFilled_),
    isAllocated_(src.isAllocated_),
    localOffset_(src.localOffset_),
    localSize_(src.localSize_),
    amap_(src.amap_),
    Amat_(NULL),
    arraysAllocated_(src.arraysAllocated_),
    val(NULL),
    bindx(NULL),
    rowLengths_(NULL),
    nnzeros_(src.nnzeros_),
    N_update_(src.N_update_),
    tmp_array_(0),
    tmp_array_len_(0),
    dtmp_array_(0),
    dtmp_array_len_(0),
    azTransformed_(src.azTransformed_)
  {
    expand_array(tmp_array_, tmp_array_len_, src.tmp_array_len_);
    expand_array(dtmp_array_, dtmp_array_len_, src.dtmp_array_len_);

    if (N_update_ > 0) {
      rowLengths_ = new int[N_update_];
      for(int i=0; i<N_update_; i++) {
        rowLengths_[i] = src.rowLengths_[i];
      }
    }

    Amat_ = AZ_matrix_create(N_update_);

    if (isAllocated_ && nnzeros_ > 0) {
      val = new double[nnzeros_+1];
      bindx = new int[nnzeros_+1];

      for(int i=0; i<nnzeros_+1; ++i) {
        val[i] = src.val[i];
        bindx[i] = src.bindx[i];
      }

      if (isFilled_) {
        AZ_set_MSR(Amat_, bindx, val, amap_->data_org, 0, NULL, AZ_LOCAL);
      }
    }
  }

  //==============================================================================
  AztecDMSR_Matrix::~AztecDMSR_Matrix()
  {
    if (arraysAllocated_) {
      delete [] val;
      val = NULL;
      delete [] bindx;
      bindx = NULL;
      arraysAllocated_ = false;
    }

    if (dtmp_array_) delete [] dtmp_array_; dtmp_array_ = NULL;

    if (N_update_ > 0) {
      delete [] rowLengths_;
      rowLengths_ = NULL;
      N_update_ = 0;
    }

    if (azTransformed_) {
      azTransformed_ = false;
    }

    AZ_matrix_destroy(&Amat_);
    Amat_ = NULL;

    isFilled_ = false;
    isAllocated_ = false;
    arraysAllocated_ = false;

    if (tmp_array_len_ > 0) {
      delete [] tmp_array_;
      tmp_array_ = 0;
      tmp_array_len_ = 0;
    }
  }

  //==============================================================================
  void AztecDMSR_Matrix::expand_array(int*& array, int& arraylen, int newlen)
  {
    if (arraylen < newlen) {
      delete [] array;
      array = new int[newlen];
      for(int i=0; i<newlen; ++i) array[i] = -999;
      arraylen = newlen;
    }
  }

  //==============================================================================
  void AztecDMSR_Matrix::expand_array(double*& array, int& arraylen, int newlen)
  {
    if (arraylen < newlen) {
      delete [] array;
      array = new double[newlen];
      for(int i=0; i<newlen; ++i) array[i] = -999.9;
      arraylen = newlen;
    }
  }

  //==============================================================================
  void AztecDMSR_Matrix::scale(double scalar)
  {
    if (scalar == 1.0) return;

    if (val != NULL) {
      for(int i=0; i<nnzeros_+1; ++i) {
        val[i] *= scalar;
      }
    }
  }

  //==============================================================================
  void AztecDMSR_Matrix::matvec(const Aztec_LSVector& x,
      Aztec_LSVector& y) const
  {
    //
    // This function forms the product y = Ax
    //

    assert(isFilled());

    int *proc_config = amap_->getProcConfig();
    //    int *idummy = 0, one=1;
    double *b = (double*)x.startPointer();
    double *c = (double*)y.startPointer();

    AZ_MSR_matvec_mult(b, c, Amat_, proc_config);
    //val,idummy,bindx,idummy,idummy,idummy,b,c,one, data_org_);

    return;
  }

  //==============================================================================
  void AztecDMSR_Matrix::put(double s)
  {
    if(isAllocated()){
      for(int i=0; i<bindx[N_update_]; i++) val[i] = s;
    }
    else {
      fei::console_out() << "AztecDMSR_Matrix::put - ERROR, can't do put until allocated"
        << FEI_ENDL;
    }
    return;
  }

  //==============================================================================
  int AztecDMSR_Matrix::rowLength(int row) const
  {
    int localRow;

    int thisRow = row;

    if(amap_->inUpdate(thisRow,localRow)){
      return(ADMSR_LOCAL_ROW_ALLOC_LEN(localRow));
    }
    else {
      fei::console_out() << "AztecDMSR_Matrix::rowLength: ERROR row " << row
        << " not in local update set." << FEI_ENDL;
      abort();
      return(-1);
    }
  }

  /**==========================================================================**/
  int AztecDMSR_Matrix::setDiagEntry(int row, double value)
  {
    int thisRow = row;
    int localRow = -1;
    if(!amap_->inUpdate(thisRow,localRow)){
      fei::console_out() << "AztecDMSR_Matrix::setDiagEntry: ERROR - row " << row
        << " not in local update set." << FEI_ENDL;
      abort(); return(-1);
    }

    val[localRow] = value;
    return(0);
  }

  /**==========================================================================**/
  double AztecDMSR_Matrix::getDiagEntry(int row) const
  {
    int thisRow = row;
    int localRow = -1;
    if(!amap_->inUpdate(thisRow,localRow)){
      fei::console_out() << "AztecDMSR_Matrix::getDiagEntry: ERROR - row " << row
        << " not in local update set." << FEI_ENDL;
      abort(); return val[0];
    }

    return(val[localRow]);
  }

  /**==========================================================================**/
  int AztecDMSR_Matrix::getOffDiagRowPointers(int row, int*& colIndices,
      double*& coefs,
      int& offDiagRowLength)
  {
    int thisRow = row;
    int localRow = -1;
    if(!amap_->inUpdate(thisRow,localRow)){
      fei::console_out() << "AztecDMSR_Matrix::getOffDiagRowPointers: ERROR - row " << row
        << " not in local update set." << FEI_ENDL;
      abort(); return(-1);
    }

    offDiagRowLength = rowLengths_[localRow];
    if (isFilled_) offDiagRowLength = ADMSR_LOCAL_ROW_ALLOC_LEN(localRow)-1;
    int start = bindx[localRow];
    colIndices = &(bindx[start]);
    coefs = &(val[start]);

    return(0);
  }

  /**==========================================================================**/
  void AztecDMSR_Matrix::getRow(int row,
      int &length,
      double *coefs,
      int *colInd) const
  {
    //
    //Note to myself:
    //getRow, putRow and sumIntoRow are incredibly ugly, and need re-written.
    //

    bool foundDiagonal = false;
    double dtmp;
    int j, localRow, itmp;

    int thisRow = row;

    if(!amap_->inUpdate(thisRow,localRow)){
      fei::console_out() << "AztecDMSR_Matrix::getRow: ERROR - row " << row
        << " not in local update set." << FEI_ENDL;
      length = 0;
      return;
    }

    int start = bindx[localRow];
    int end = bindx[localRow+1]-1;

    j = 0;
    for(int i=start; i<=end; i++){
      coefs[j] = val[i];
      colInd[j] = amap_->getTransformedEqn(bindx[i]);

      if (colInd[j]==row) {
        //we're at the diagonal element, so put it in.
        dtmp = coefs[j];
        itmp = colInd[j];
        coefs[j] = val[localRow];
        colInd[j] = row;
        j++;
        coefs[j] = dtmp;
        colInd[j] = itmp;
        foundDiagonal = true;
      }
      j++;
    }

    if(!foundDiagonal){ // still need to put in the diagonal
      coefs[j] = val[localRow];
      colInd[j] = row;
    }

    length = j+1;
    return;
  }

  /**==========================================================================**/
  int AztecDMSR_Matrix::putRow(int row, int len, const double *coefs,
      const int *colInd)
  {
    //
    //Note to myself:
    //getRow, putRow and sumIntoRow are incredibly ugly, and need re-written.
    //

    int j, localRow, globalColIndex;

    int thisRow = row;

    if (!amap_->inUpdate(thisRow,localRow)){
      fei::console_out() << "AztecDMSR_Matrix::putRow: ERROR row " << row
        << " not in local update set." << FEI_ENDL;
      return(-1);
    }

    int offDiagRowLen = rowLengths_[localRow];
    int rowAllocLen = ADMSR_LOCAL_ROW_ALLOC_LEN(localRow);
    if (isFilled_) offDiagRowLen = rowAllocLen - 1;
    int offDiagRowAllocLen = rowAllocLen - 1;

    if (len > rowAllocLen) {
      fei::console_out() << "AztecDMSR_Matrix::putRow. too many coefs, row " << row << FEI_ENDL;
      return(-1);
    }

    int jLimit = bindx[localRow+1] - 1;
    int jStart = bindx[localRow];
    int* colInds = &(bindx[jStart]);
    int colStart = *colInds;
    double* rowCoefs = &(val[jStart]);

    for(int i=0; i<len; i++){
      if (colInd[i] == row){  //it's on the diagonal
        val[localRow] = coefs[i];
        continue;
      }

      //look along the row until we find the right col index
      //or an empty spot
      j=jStart;
      if (isFilled()){
        int colIndex = colStart;
        int col = colInd[i];
        while (j <= jLimit) {
          globalColIndex = amap_->getTransformedEqn(colIndex);

          if (globalColIndex == col) break;

          colIndex = bindx[++j];
        }

        //now put the coefficient in if we haven't gone too far
        //along the row
        if (j <= jLimit) {
          val[j] = coefs[i];
        }
        else {
          fei::console_out() << "AztecDMSR_Matrix::putRow: ERROR didn't "
            << " find col. index " << colInd[i] << " in "
            << "row " << row << FEI_ENDL;
          return(-1);
        }
      }
      else { // !isFilled()

        //first, look for colInd[i] in the row
        int col = colInd[i];
        int insertPoint = -1;
        int index = fei::binarySearch<int>(col, colInds, offDiagRowLen, insertPoint);

        if (index >= offDiagRowAllocLen){ //bad index
          fei::console_out() << "AztecDMSR_Matrix::putRow, ERROR: "
            << "row " << row << ", colInd["<<i<<"] " << colInd[i]
            << ", index = " << index << FEI_ENDL;
          return(-1);
        }

        if (index >= 0) {
          rowCoefs[index] = coefs[i];
        }
        else {
          int tmp = offDiagRowLen;
          int err = insert(col, insertPoint, colInds,
              tmp, offDiagRowAllocLen);
          err += insert(coefs[i], insertPoint, rowCoefs,
              offDiagRowLen, offDiagRowAllocLen);
          if (err != 0) {
            fei::console_out() << "AztecDMSR_Matrix::putRow ERROR: failed to add "
              << "value for index " << col << " to row " << row << FEI_ENDL;
            return(-1);
          }
          rowLengths_[localRow]++;
        }
      }
    }

    return(0);
  }

  /**==========================================================================**/
  int AztecDMSR_Matrix::sumIntoRow(int numRows, const int* rows,
      int numCols, const int *colInd,
      const double* const* coefs)
  {
    if (numRows == 0 || numCols == 0) return(0);

    if (!isFilled_) {
      int err = 0;
      for(int i=0; i<numRows; ++i) {
        err = sumIntoRow(rows[i], numCols, coefs[i], colInd);
        if (err != 0) return(err);
      }

      return(0);
    }

    //Now for the harder (but more important) case where isFilled_ == true.

    //first compute max-row-length:
    int maxRowLen = 0;
    for(int i=0; i<numRows; ++i) {
      int row = rows[i];
      int localRow;
      if (!amap_->inUpdate(row, localRow)) {
        fei::console_out() << "AztecDMSR_Matrix::sumIntoRow: ERROR row " << row
          << " not in local update set [" << amap_->getUpdate()[0] << " ... "
          << amap_->getUpdate()[N_update_-1] << "]." << FEI_ENDL;
        return(-1);
      }

      int rowlen = bindx[localRow+1]-bindx[localRow];
      if (maxRowLen < rowlen) maxRowLen = rowlen;
    }

    if (maxRowLen+2*numCols > tmp_array_len_) {
      expand_array(tmp_array_, tmp_array_len_, maxRowLen+2*numCols);
    }

    int* incols = &tmp_array_[maxRowLen];
    int* indirect = incols+numCols;

    for(int jj=0; jj<numCols; ++jj) {
      incols[jj] = colInd[jj];
      indirect[jj] = jj;
    }

    fei::insertion_sort_with_companions<int>(numCols, incols, indirect);

    int row, localRow;

    for(int i=0; i<numRows; ++i) {
      row = rows[i];
      if (!amap_->inUpdate(row, localRow)) {
        fei::console_out() << "AztecDMSR_Matrix::sumIntoRow: ERROR row " << row
          << " not in local update set [" << amap_->getUpdate()[0] << " ... "
          << amap_->getUpdate()[N_update_-1] << "]." << FEI_ENDL;
        return(-1);
      }

      int jStart = bindx[localRow];
      double* rowCoefs = val+jStart;
      int* rowColInds = bindx+jStart;
      int rowLen= bindx[localRow+1]-jStart;

      for(int jj=0; jj<rowLen; ++jj) {
        tmp_array_[jj] = amap_->getTransformedEqn(rowColInds[jj]);
      }

      const double* coefs_i = coefs[i];

      int inoffset = 0;
      int incol = incols[inoffset];
      while (incol == row) {
        val[localRow] += coefs_i[indirect[inoffset++]];
        if (inoffset >= numCols) break;
        incol = incols[inoffset];
      }

      if (inoffset >= numCols) continue;

      //rowOffset is the offset into the row at which incol appears.

      int rowOffset = fei::binarySearch<int>(incol, tmp_array_, rowLen);
      if (rowOffset < 0) {
        fei::console_out() << "AztecDMSR_Matrix::sumIntoRow, ERROR: "
          << "row " << row << ", col not found: "
          << incol << FEI_ENDL;
        return(-1);
      }

      rowCoefs[rowOffset++] += coefs_i[indirect[inoffset++]];

      //check whether incols has a repeated column-index
      if (inoffset>0 && incols[inoffset] == incols[inoffset-1]) --rowOffset;

      while(inoffset < numCols) {
        incol = incols[inoffset];

        if (incol == row) {
          val[localRow] += coefs_i[indirect[inoffset++]];
          continue;
        }

        while(tmp_array_[rowOffset] != incol) {
          ++rowOffset;
          if (rowOffset >= rowLen) {
            fei::console_out() << "AztecDMSR_Matrix::sumIntoRow, ERROR, col "
              << incol << " not found in row " << row << FEI_ENDL;
            return(-1);
          }
        }

        rowCoefs[rowOffset++] += coefs_i[indirect[inoffset++]];
        if (inoffset>0 && incols[inoffset] == incols[inoffset-1]) --rowOffset;
      }
    }

    return(0);
  }

  /**==========================================================================**/
  int AztecDMSR_Matrix::sumIntoRow(int row, int len, const double *coefs,
      const int *colInd)
  {
    //
    //Note to myself:
    //getRow, putRow and sumIntoRow are incredibly ugly, and need re-written.
    //

    int localRow, thisRow = row ;

    if (!amap_->inUpdate(thisRow,localRow)) {
      fei::console_out() << "AztecDMSR_Matrix::sumIntoRow: ERROR row " << row
        << " not in local update set." << FEI_ENDL;
      return(-1);
    }

    int jLimit = bindx[localRow+1] - 1;
    int jStart = bindx[localRow];
    int jLen = jLimit-jStart+1;
    int* colInds = &(bindx[jStart]);
    double* rowCoefs = &(val[jStart]);

    if (isFilled_) {
      if (jLen+len > tmp_array_len_) {
        expand_array(tmp_array_, tmp_array_len_, jLen+len);
        expand_array(dtmp_array_, dtmp_array_len_, len);
      }
      for(int jj=0; jj<jLen; ++jj) {
        tmp_array_[jj] = amap_->getTransformedEqn(colInds[jj]);
      }

      int* incols = &tmp_array_[jLen];
      int doffs = 0;
      for(int jj=0; jj<len; ++jj) {
        int col = colInd[jj];
        if (col == row) {
          val[localRow] += coefs[jj];
        }
        else {
          incols[doffs] = col;
          dtmp_array_[doffs++] = coefs[jj];
        }
      }
      fei::insertion_sort_with_companions<double>(doffs, incols, dtmp_array_);

      int ioffset = 0;
      int offset = fei::binarySearch<int>(incols[ioffset], tmp_array_, jLen);
      if (offset < 0) {
        fei::console_out() << "AztecDMSR_Matrix::sumIntoRow, ERROR: "
          << "row " << row << ", col not found: "
          << colInd[ioffset] << FEI_ENDL;
        return(-1);
      }

      rowCoefs[offset++] += dtmp_array_[ioffset++];
      if (incols[ioffset] == tmp_array_[offset-1]) --offset;

      while(ioffset < doffs) {
        int incol = incols[ioffset];

        while(tmp_array_[offset] != incol) {
          ++offset;
          if (offset >= jLen) {
            fei::console_out() << "AztecDMSR_Matrix::sumIntoRow, ERROR, col "
              << incols[ioffset] << " not found in row " << row << FEI_ENDL;
            return(-1);
          }
        }
        rowCoefs[offset++] += dtmp_array_[ioffset++];
        if (incols[ioffset] == tmp_array_[offset-1]) --offset;
      }

      return(0);
    }

    //if we get to here, then we know that isFilled_ is false...

    int rowAllocLen = ADMSR_LOCAL_ROW_ALLOC_LEN(localRow);
    int offDiagRowLen = isFilled_ ? rowAllocLen - 1 : rowLengths_[localRow];
    int offDiagRowAllocLen = rowAllocLen - 1;

    for(int i=0; i<len; i++){
      if (colInd[i] == row){  //it's on the diagonal
        val[localRow] += coefs[i];
        continue;
      }

      //find the right col index in the row, or an empty spot
      int col = colInd[i];
      int insertPoint = -1;
      int index = fei::binarySearch<int>(col, colInds, offDiagRowLen, insertPoint);

      if (index >= 0) {
        rowCoefs[index] += coefs[i];
      }
      else {
        int tmp = offDiagRowLen;
        int err = insert(col, insertPoint, colInds,
            tmp, offDiagRowAllocLen);
        err += insert(coefs[i], insertPoint, rowCoefs,
            offDiagRowLen, offDiagRowAllocLen);
        if (err != 0) {
          fei::console_out() << "AztecDMSR_Matrix::sumIntoRow ERROR: failed to add "
            << "value for index " << col << " to row " << row << FEI_ENDL;
          return(-1);
        }
        rowLengths_[localRow]++;
      }
    }

    return(0);
  }

  //==============================================================================
  int AztecDMSR_Matrix::addScaledMatrix(double scalar,
      const AztecDMSR_Matrix& source)
  {
    if (N_update_ != source.N_update_ ||
        nnzeros_ != source.nnzeros_ ||
        isFilled_ != source.isFilled_) {
      fei::console_out() << "AztecDMSR_Matrix::addScaledMatrix ERROR, not compatible"
        << FEI_ENDL;
      return(-1);
    }

    const double* src_val = source.val;
    int i;
    for(i=0; i<N_update_; ++i) {
      val[i] += scalar*src_val[i];
    }

    //val[N_update_] is not used.

    if (scalar == 1.0) {
      for(i=N_update_+1; i<nnzeros_+1; ++i) {
        val[i] += src_val[i];
      }
    }
    else {
      for(i=N_update_+1; i<nnzeros_+1; ++i) {
        val[i] += scalar*src_val[i];
      }
    }

    return(0);
  }

  //==============================================================================
  int AztecDMSR_Matrix::insert(int item, int offset, int* list,
      int& len, int allocLen)
  {
    if (len >= allocLen) return(-1);

    for(int i=len; i>offset; i--) list[i] = list[i-1];

    list[offset] = item;
    len++;

    return(0);
  }

  //==============================================================================
  int AztecDMSR_Matrix::insert(double item, int offset, double* list,
      int& len, int allocLen)
  {
    if (len >= allocLen) return(-1);

    for(int i=len; i>offset; i--) list[i] = list[i-1];

    list[offset] = item;
    len++;

    return(0);
  }

  //==============================================================================
  void AztecDMSR_Matrix::getDiagonal(Aztec_LSVector& diagVector) const {

    /** AztecDMSR_Matrix::getDiagonal --- form a Vector of the diagonals of
      the matrix **/

    // have each processor form its piece of the diagonal
    double *pdv = (double*)diagVector.startPointer();

    for(int i=0; i<N_update_; i++){
      pdv[i] = val[i];
    }
  }

  /**==========================================================================**/
  void AztecDMSR_Matrix::allocate(int *rowLengths)
  {
    //
    //We assume that 'rowLengths' is of length N_update_.
    //
    //rowLengths contains the length of each local row, *NOT* including the
    //coefficient on the diagonal.
    //
    int i;

    //first, count how many non-zeros there are in the local submatrix

    nnzeros_ = 0;
    for(i=0; i<N_update_; i++){
      if (rowLengths[i] < 0) {
        messageAbort("allocate: negative row length");
      }
      nnzeros_ += rowLengths[i] + 1;
    }

    if (bindx != NULL) {
      delete [] bindx;
    }
    bindx = new int[nnzeros_+1];

    if (val != NULL) {
      delete [] val;
    }
    val = new double[nnzeros_+1];

    arraysAllocated_ = true;

    for(i=0; i<nnzeros_+1; i++){
      val[i] = 0.0;
      bindx[i] = -1;
    }

    bindx[0] = N_update_+1;

    for(i=0; i<N_update_; i++){
      bindx[i+1] = bindx[i] + rowLengths[i];
      if (bindx[i+1] < 0) {
        messageAbort("allocate: bindx row length negative.");
      }
    }

    //val[N_update_] not used by aztec but we'll initialize it anyway...
    val[N_update_] = 0.0;

    AZ_set_MSR(Amat_, bindx, val,amap_->data_org, 0, NULL, AZ_LOCAL);

    setAllocated(true);
    return;
  }

  /**=========================================================================**/
  void AztecDMSR_Matrix::allocate(int *rowLengths,
      const int* const* colIndices)
  {
    allocate(rowLengths);
    int col;

    int offset = N_update_+1;
    for(int i=0; i<N_update_; ++i) {
      const int* row_cols = colIndices[i];
      int row_len = rowLengths[i];
      rowLengths_[i] = row_len-1;

      int prev_col = -999;
      int coffset = 0;
      for(int j=0; j<row_len; ++j) {
        col = row_cols[coffset++];

        if (col == localOffset_+i) {
          col = row_cols[coffset++];
        }

        if (col <= prev_col) {
          messageAbort("allocate: column-indices not sorted.");
        }

        prev_col = col;

        bindx[offset++] = col;
      }
    }
  }

  /**=========================================================================**/
  double AztecDMSR_Matrix::rowMax(int row) const {
    int localRow;
    double max = 0.0;

    if(!amap_->inUpdate(row,localRow)){
      fei::console_out() << "AztecDMSR_Matrix::rowMax: ERROR row " << row
        << " not in local update set." << FEI_ENDL;
      return(-1.0);
    }

    max = fabs(val[localRow]);

    for(int i=bindx[localRow]; i<bindx[localRow+1]; i++)
      if(fabs(val[i])>max)max = fabs(val[i]);

    return(max);
  }

  /**==========================================================================**/
  void AztecDMSR_Matrix::fillComplete() {
    /*
       This is where we call the Aztec function AZ_transform, which calculates
       communication parameters and re-orders the equations for use as a
       global distributed matrix.
       */
    if (isFilled_ || azTransformed_) {
      isFilled_ = true;
      azTransformed_ = true;
      return;
    }

    int *proc_config = amap_->getProcConfig();
    int *dummy = 0;

    //before we turn Aztec loose on the matrix, lets do a quick check on the
    //indices to try to make sure none of them are garbage...
    int globalSize = amap_->globalSize();
    for(int i=N_update_+1; i<nnzeros_+1; i++) {
      if (bindx[i] < 0 || bindx[i] >= globalSize) {
        fei::console_out() << "AztecDMSR_Matrix: ERROR, bindx["<<i<<"]: " << bindx[i]
          << ", globalSize: " << globalSize << FEI_ENDL;
#ifndef FEI_SER
        MPI_Comm thisComm = amap_->getCommunicator();
        MPI_Abort(thisComm, -1);
#endif
      }
    }

    AZ_transform(proc_config, &amap_->external, bindx, val,
        amap_->getUpdate(), &amap_->update_index,
        &amap_->extern_index, &amap_->data_org, N_update_,
        dummy, dummy, dummy, &dummy, AZ_MSR_MATRIX);

    //AZ_transform allocates these arrays:
    //  amap_->external
    //  amap_->update_index
    //  amap_->extern_index
    //  amap_->data_org
    //
    //On return from AZ_transform, the array update_index contains a mapping
    //to the local re-ordering of the indices of the update array. Now we will fill
    //the orderingUpdate array with the reverse of that mapping. i.e., a record
    //of how to get back to the original ordering of the update indices.

    AZ_set_MSR(Amat_, bindx, val, amap_->data_org, 0, NULL, AZ_LOCAL);

    amap_->orderingUpdate.resize(N_update_);
    for(int ii=0; ii<N_update_; ii++) {
      amap_->orderingUpdate[amap_->update_index[ii]] = ii;
    }

    amap_->az_transformed = true;
    azTransformed_ = true;

    setFilled(true);
    return;
  }

  //==============================================================================
  void AztecDMSR_Matrix::copyStructure(AztecDMSR_Matrix& source)
  {
    //
    //This function copies the structure (essentially just the bindx and
    //rowLengths_ arrays) and other relevant variables from the 'source' matrix.
    //The result is that 'this' matrix is laid out the same as 'source'.
    //
    nnzeros_ = source.nnzeros_;

    if (arraysAllocated_) {
      delete [] val;
      delete [] bindx;
      arraysAllocated_ = false;
    }

    val = new double[nnzeros_+1];
    bindx = new int[nnzeros_+1];

    int i;
    for(i=0; i<nnzeros_+1; i++) {
      val[i] = 0.0;
      bindx[i] = source.bindx[i];
    }

    for(i=0; i<N_update_; ++i) rowLengths_[i] = source.rowLengths_[i];

    amap_ = source.amap_;

    AZ_set_MSR(Amat_, bindx, val, amap_->data_org, 0, NULL, AZ_LOCAL);

    isFilled_ = source.isFilled_;
    azTransformed_ = source.azTransformed_;

    arraysAllocated_ = true;
    setAllocated(true);
  }

  //=============================================================================
  bool AztecDMSR_Matrix::readFromFile(const char *filename)
  {
    /*
       This function reads the matrix data from a matrix-market format data file,
       which looks like:

       n
       i j val
       .
       .

       Important note: we are going to assume that the index-base of the indices in
       the file are 0.
       */
    int i, j, dummy;
    double value;
    char line[128];

    FILE *mfp = fopen(filename,"r");

    if(!mfp){
      fei::console_out() << "AztecDMSR_Matrix::readFromFile - couldn't open matrix file."
        << FEI_ENDL;
      return(false);
    }

    if (strstr(filename, ".mtx") == NULL) {
      fei::console_out() << "AztecDMSR_Matrix::readFromFile: filename doesn't contain "
        << "'.mtx'. File should be a MatrixMarket file." << FEI_ENDL;
      return(false);
    }

    do {
      fgets(line,128,mfp);
    } while(strchr(line,'%'));

    while(!feof(mfp)){
      do {
        fgets(line,128,mfp);
      } while(strchr(line,'%'));
      if(feof(mfp)){
        fclose(mfp);
        return(true);
      }
      sscanf(line,"%d %d %le",&i,&j,&value);

      if(amap_->inUpdate(i, dummy)) {
        if (putRow(i, 1, &value, &j) != 0) return(false);
      }
    }

    fclose(mfp);
    return(true);
  }

  /**=========================================================================**/
  bool AztecDMSR_Matrix::writeToFile(const char *fileName) const
  {
    /* Write the matrix into the file "fileName", using the format:
       n
       i j val
       .
       .
       */

    int numProcs = amap_->getProcConfig()[AZ_N_procs];
    int thisProc = amap_->getProcConfig()[AZ_node];
    int masterRank = 0;


    int localNNZ = nnzeros_;
    int globalNNZ = localNNZ;
#ifndef FEI_SER
    MPI_Comm thisComm = amap_->getCommunicator();
    MPI_Allreduce(&localNNZ, &globalNNZ, 1, MPI_INT, MPI_SUM, thisComm);
#endif

    for(int p=0; p<numProcs; p++){

      //A barrier inside the loop so each processor waits its turn.
#ifndef FEI_SER
      MPI_Barrier(thisComm);
#endif

      if (p == thisProc) {
        FILE *file = NULL;

        if (masterRank == thisProc) {
          //This is the master processor, open a new file.
          file = fopen(fileName,"w");

          //Write the matrix dimensions n and n (rows==cols) into the file,
          //along with the global number-of-nonzeros globalNNZ.

          int n = amap_->globalSize();
          fprintf(file,"%d %d %d\n",n, n, globalNNZ);
        }
        else {
          //This is a non-master node, open the file for appending to.
          file = fopen(fileName,"a");
        }

        //Now loop over the local portion of the matrix.
        for(int i=0; i<N_update_; i++){
          int row = localOffset_+i;

          int localRow = -1;
          if (!amap_->inUpdate(row, localRow)) return(false);

          int offDiagRowLen = ADMSR_LOCAL_ROW_LEN(localRow) - 1;
          if (isFilled_) offDiagRowLen = ADMSR_LOCAL_ROW_ALLOC_LEN(localRow) - 1;
          int* colInds = &(bindx[bindx[localRow]]);
          double* coefs = &(val[bindx[localRow]]);

          bool wroteDiagonal = false;
          for(int j=0; j<offDiagRowLen; j++) {
            int col = colInds[j];
            int globalCol = col;
            if (isFilled()) {
              globalCol = amap_->getTransformedEqn(col);
            }

            if (globalCol >= row && !wroteDiagonal) {
              fprintf(file,"%d %d %20.13e\n", row, row,
                  val[localRow]);
              wroteDiagonal = true;
            }

            fprintf(file,"%d %d %20.13e\n", row, globalCol, coefs[j]);
          }
          if (!wroteDiagonal) {
            fprintf(file,"%d %d %20.13e\n", row, row,
                val[localRow]);
            wroteDiagonal = true;
          }
        }
        fclose(file);
      }
    }

    return(true);
  }

  //==============================================================================
  void AztecDMSR_Matrix::messageAbort(const char* mesg) {
    fei::console_out() << "AztecDMSR_Matrix: ERROR: " << mesg << " Aborting." << FEI_ENDL;
    abort();
  }

}//namespace fei_trilinos

#endif
//HAVE_FEI_AZTECOO

