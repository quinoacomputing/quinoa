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


#include <fei_FillableMat.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_CSVec.hpp>

namespace fei {

//-----------------------------------------------------------------
FillableMat::FillableMat()
 : matdata_(),
   vecpool_()
{
}

//-----------------------------------------------------------------
FillableMat::FillableMat(EqnBuffer& eqnbuf)
 : matdata_(),
   vecpool_()
{
  std::vector<int>& eqnNums = eqnbuf.eqnNumbers();
  int numEqns = eqnNums.size();
  std::vector<fei::CSVec*>& eqns = eqnbuf.eqns();

  for(int i=0; i<numEqns; ++i) {
    int row = eqnNums[i];
    fei::CSVec* row_vec = eqns[i];
    int rowlen = row_vec->size();
    int* indices = &(row_vec->indices()[0]);
    double* coefs = &(row_vec->coefs()[0]);

    for(int j=0; j<rowlen; ++j) {
      putCoef(row, indices[j], coefs[j]);
    }
  }
}

//-----------------------------------------------------------------
FillableMat::~FillableMat()
{
  feipoolmat::iterator
    iter = matdata_.begin(), iter_end = matdata_.end();
  for(; iter!=iter_end; ++iter) {
    vecpool_.destroy(iter->second);
    vecpool_.deallocate(iter->second, 1);
  }
}

//-----------------------------------------------------------------
FillableMat&
FillableMat::operator=(const FillableMat& src)
{
  clear();

  FillableMat::const_iterator
    s_iter = src.begin(),
    s_end = src.end();

  for(; s_iter != s_end; ++s_iter) {
    int row = s_iter->first;
    const FillableVec* srow = s_iter->second;

    FillableVec::const_iterator
      r_iter = srow->begin(),
      r_end = srow->end();

    for(; r_iter != r_end; ++r_iter) {
      int col = r_iter->first;
      double coef = r_iter->second;

      putCoef(row, col, coef);
    }
  }

  return *this;
}

//-----------------------------------------------------------------
void
FillableMat::setValues(double value)
{
  feipoolmat::iterator
    iter = matdata_.begin(), iter_end = matdata_.end();

  for(; iter != iter_end; ++iter) {
    iter->second->setValues(value);
  }
}

//-----------------------------------------------------------------
void
FillableMat::createPosition(int row, int col)
{
  sumInCoef(row, col, 0.0);
}

//-----------------------------------------------------------------
FillableMat::feipoolmat::iterator
insert_row(FillableMat::feipoolmat& matdata,
           FillableMat::feipoolmat::iterator iter,
           int row,
           fei_Pool_alloc<FillableVec>& vecpool)
{
  static FillableVec dummy;

  FillableVec* vptr = vecpool.allocate(1);
  vecpool.construct(vptr, dummy);

  return matdata.insert(iter, std::make_pair(row, vptr));
}

//-----------------------------------------------------------------
void
FillableMat::sumInCoef(int row, int col, double coef)
{
  FillableVec* rowvec = create_or_getRow(row);

  rowvec->addEntry(col, coef);
}

//-----------------------------------------------------------------
void
FillableMat::putCoef(int row, int col, double coef)
{
  FillableVec* rowvec = create_or_getRow(row);

  rowvec->putEntry(col, coef);
}

//-----------------------------------------------------------------
void
FillableMat::sumInRow(int row, const int* cols, const double* coefs,
                      unsigned len)
{
  FillableVec* rowvec = create_or_getRow(row);

  for(unsigned i=0; i<len; ++i) {
    rowvec->addEntry(cols[i], coefs[i]);
  }
}

//-----------------------------------------------------------------
void
FillableMat::putRow(int row, const int* cols, const double* coefs,
                    unsigned len)
{
  FillableVec* rowvec = create_or_getRow(row);

  for(unsigned i=0; i<len; ++i) {
    rowvec->putEntry(cols[i], coefs[i]);
  }
}

//-----------------------------------------------------------------
unsigned
FillableMat::getNumRows() const
{
  return matdata_.size();
}

//-----------------------------------------------------------------
bool
FillableMat::hasRow(int row) const
{
  feipoolmat::const_iterator iter = matdata_.find(row);
  return iter != matdata_.end();
}

//-----------------------------------------------------------------
const FillableVec*
FillableMat::getRow(int row) const
{
  feipoolmat::const_iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    throw std::runtime_error("fei::FillableMat: row not found.");
  }

  return iter->second;
}

//-----------------------------------------------------------------
FillableVec*
FillableMat::create_or_getRow(int row)
{
  feipoolmat::iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    iter = insert_row(matdata_, iter, row, vecpool_);
  }

  return iter->second;
}

//-----------------------------------------------------------------
void
FillableMat::clear()
{
  feipoolmat::iterator
    iter = matdata_.begin(), iter_end = matdata_.end();
  for(; iter!=iter_end; ++iter) {
    vecpool_.destroy(iter->second);
    vecpool_.deallocate(iter->second, 1);
  }

  matdata_.clear();
}

//-----------------------------------------------------------------
bool
FillableMat::operator==(const FillableMat& rhs) const
{
  if (getNumRows() != rhs.getNumRows()) return false;

  FillableMat::const_iterator
    this_it = begin(),
    this_end = end();

  FillableMat::const_iterator
    rhs_it = rhs.begin(),
    rhs_end = rhs.end();

  for(; this_it != this_end; ++this_it, ++rhs_it) {
    int this_row = this_it->first;
    int rhs_row = rhs_it->first;
    if (this_row != rhs_row) return false;

    const FillableVec* this_row_vec = this_it->second;
    const FillableVec* rhs_row_vec = rhs_it->second;

    if (*this_row_vec != *rhs_row_vec) return false;
  }

  return true;
}

//-----------------------------------------------------------------
bool
FillableMat::operator!=(const FillableMat& rhs) const
{
  return !(*this == rhs);
}

//-----------------------------------------------------------------
void print(std::ostream& os, const FillableMat& mat)
{
  FillableMat::const_iterator
    irow = mat.begin(), irowend = mat.end();
  for(; irow!=irowend; ++irow) {
    int row = irow->first;
    const FillableVec* vec = irow->second;
    os << "row " << row << ": ";
    FillableVec::const_iterator
      ivec = vec->begin(), ivecend = vec->end();
    for(; ivec!=ivecend; ++ivec) {
      os << "("<<ivec->first<<","<<ivec->second<<") ";
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------
int count_nnz(const FillableMat& mat)
{
  int nnz = 0;

  FillableMat::const_iterator
    r_iter = mat.begin(),
    r_end = mat.end();

  for(; r_iter != r_end; ++r_iter) {
    FillableVec* row = r_iter->second;
    nnz += row->size();
  }

  return nnz;
}

//-----------------------------------------------------------------
void get_row_numbers(const FillableMat& mat, std::vector<int>& rows)
{
  rows.resize(mat.getNumRows());

  FillableMat::const_iterator
    m_iter = mat.begin(),
    m_end = mat.end();

  size_t offset = 0;
  for(; m_iter!=m_end; ++m_iter) {
    rows[offset++] = m_iter->first;
  }
}

}//namespace fei

