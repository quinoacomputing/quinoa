// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef COO_MATRIX_PARTITIONED_VIEW_CLASS_DEF_H
#define COO_MATRIX_PARTITIONED_VIEW_CLASS_DEF_H

#include <sstream>
#include <algorithm>
#include <functional>

#include "AbstractLinAlgPack_COOMatrixPartitionedViewClassDecl.hpp"

namespace AbstractLinAlgPack {

// /////////////////////////////////////////////////////////////////////////////
// Template function definitions for COOMatrixPartitionedView

template <class T_Indice, class T_Value>
void COOMatrixPartitionedView<T_Indice,T_Value>::create_view(
        size_type				rows
      , size_type				cols
      , size_type				nz		
      , value_type			val[]
      , const indice_type		ivect[]
      , const indice_type		jvect[]
      , const size_type		inv_row_perm[]
      , const size_type		inv_col_perm[]
      , const size_type		num_row_part
      , const size_type		row_part[]
      , const size_type		num_col_part
      , const size_type		col_part[]
      , const EPartitionOrder	partition_order )
{
  try {

  // 1) Check some preconditins before we start

  // Check the sparsisty density of COO matrix
  if(nz > rows * cols)
    throw std::out_of_range(
      "COOMatrixPartitionedView<...>::create_view() : Input error, "
      "nz can not be greater than rows * cols");
  
  // Check row and column partition information

  if(num_row_part > rows || num_col_part > rows)
    throw std::out_of_range(
      "COOMatrixPartitionedView<...>::create_view() : Input error, "
      "num_rows_part and num_col_part can not be greater than"
      " rows and cols respectivly");

  if(row_part[num_row_part] > rows + 1)
    throw std::out_of_range(
      "COOMatrixPartitionedView<...>::create_view() : Input error, "
      "row_part[num_row_part] can not be greater than rows");

  if(col_part[num_col_part] > cols + 1)
    throw std::out_of_range(
      "COOMatrixPartitionedView<...>::create_view() : Input error, "
      "col_part[num_col_part] can not be greater than cols");

  {for(size_type i = 1; i < num_row_part + 1; ++i)
    if(row_part[i-1] >= row_part[i])
      throw std::domain_error(
        "COOMatrixPartitionedView<...>::create_view() : Input error, "
        "row_part[i-1] < row_part[i] does not hold");}

  {for(size_type i = 1; i < num_col_part + 1; ++i)
    if(col_part[i-1] >= col_part[i])
      throw std::domain_error(
        "COOMatrixPartitionedView<...>::create_view() : Input error, "
        "col_part[i-1] < col_part[i] does not hold");}
                
  // Get references to referenced quantities
  std::vector<size_type>
    &_row_part		= ref_row_part_.obj(),
    &_col_part		= ref_col_part_.obj(),
    &_part_start	= ref_part_start_.obj();
  ele_type
    &_ele			= ref_ele_.obj();

  // 2) Initialize storage for data members and copy data
  num_row_part_		= num_row_part;
  num_col_part_		= num_col_part;
  _row_part.resize(num_row_part_ + 1);
  std::copy(row_part, row_part+ num_row_part_+1, _row_part.begin());
  _col_part.resize(num_col_part_ + 1);
  std::copy(col_part, col_part+ num_col_part_+1, _col_part.begin());
  partition_order_	= partition_order;
  _ele.resize(nz,element_type());	// hack to get around compiler error.
  _part_start.resize(num_row_part_ * num_col_part_+1);
  _part_start.assign(_part_start.size(),0);	// set to 0

  // 3) Count the number of nonzero elements in each overall partition
  {
    // Use the storage locations _part_start[1] ... _part_start[num_part]
    // to count the number of nonzero elements in each overall partition.
    // In particular num_in_part[overall_p - 1] will give the number
    // of nonzero elements in the overall partition number overall_p.
    //
    // _part_start = { 0, nz1, nz2,...,nzn } 
    //
    size_type *num_in_part = &_part_start[0] + 1; 

    // Loop (time = O(nz), space = O(1))

    // read iterators
    const indice_type	*ivect_itr		= ivect,
              *ivect_itr_end	= ivect + nz,
              *jvect_itr		= jvect;

    for(;ivect_itr != ivect_itr_end; ++ivect_itr, ++jvect_itr) {
      // Get row and column indices in the non-permuted matrix
      indice_type		i_org	= *ivect_itr,
              j_org	= *jvect_itr;
      // assert that they are in range
      if(i_org < 1 || i_org > rows || j_org < 1 || j_org > cols) {
        std::ostringstream omsg;
        omsg	<<	"COOMatrixPartitionedView<...>::create_view() : "
              " Error, element k = " << ivect_itr - ivect
            <<	" in the non-permuted matrix"
              " is out of range with rows = " << rows
            <<	", cols = " << cols << ", i = " << i_org
            <<	", and j = " << j_org;
        throw std::out_of_range(omsg.str());
      }
      // Get row and column indices in the permuted matrix
      indice_type		i = inv_row_perm[i_org - 1],
              j = inv_col_perm[j_org - 1];
      // assert that they are in range
      if(i < 1 || i > rows || j < 1 || j > cols) {
        std::ostringstream omsg;
        omsg	<<	"COOMatrixPartitionedView<...>::create_view() : "
              " Error, element k = " << ivect_itr - ivect
            <<	" in the permuted matrix"
              " is out of range with rows = " << rows
            <<	", cols = " << cols << ", i = " << i
            <<	", and j = " << j;
        throw std::out_of_range(omsg.str());
      }
      // get the overall partition number
      size_type overall_p = overall_p_from_ij(_row_part,_col_part,i,j);
      // Increment the number of nonzero elements.
      num_in_part[overall_p - 1]++;
    }
  }

  // 4) Set part_start[ovarall_p] equal to the start in ele
  // for the nonzero elements in that partition.
  //
  // _part_start = { start_1 = 0, start_2,..., start_n, ###}
  //
  {for(size_type i = 2; i < num_row_part_ * num_col_part_; ++i)
    _part_start[i] += _part_start[i-1];}

  // 5) Shift the elements over
  //
  // _part_start = { 0, start_1 = 0, start_2,..., start_n}
  //
  {for(size_type i = num_row_part_ * num_col_part_; i > 0; --i)
    _part_start[i] = _part_start[i-1];}

  // 6) Add the nonzero elements to each partition.  When we
  // are done we should have _part_start initialized properly.
  //
  // part_start = { start_1 = 0, start_2,..., start_n, total_nz }
  // 
  {
    // next_ele_insert[overall_p - 1] is the possition in ele
    // for the next element to ensert
    size_type	*next_ele_insert = &_part_start[0] + 1;

    // Loop (time = O(nz), space = O(1))

    // read iterators
    value_type			*val_itr		= val,
              *val_itr_end	= val + nz;
    const indice_type	*ivect_itr		= ivect,
              *jvect_itr		= jvect;

    for(;val_itr != val_itr_end; ++val_itr, ++ivect_itr, ++jvect_itr) {
      // Get row and column indices in the permuted matrix
      indice_type		i = inv_row_perm[*ivect_itr - 1],
              j = inv_col_perm[*jvect_itr - 1];
      // get the overall partition number
      size_type overall_p = overall_p_from_ij(_row_part,_col_part,i,j);
      // Add the element to the partition
      _ele[next_ele_insert[overall_p - 1]++].initialize(val_itr,i,j);
    }
  }

  } // end try
  catch(...) {
    free();
    throw;	// rethrow the exception out of here
  }
}

template <class T_Indice, class T_Value>
void COOMatrixPartitionedView<T_Indice,T_Value>::bind(const COOMatrixPartitionedView& coo_view)
{
  num_row_part_		= coo_view.num_row_part_;
  num_col_part_		= coo_view.num_col_part_;
  ref_row_part_		= coo_view.ref_row_part_;
  ref_col_part_		= coo_view.ref_col_part_;
  partition_order_	= coo_view.partition_order_;
  ref_ele_			= coo_view.ref_ele_;
  ref_part_start_		= coo_view.ref_part_start_;
}

template <class T_Indice, class T_Value>
void COOMatrixPartitionedView<T_Indice,T_Value>::free()
{
  // Reinitialize to uninitizlied
  num_row_part_ = num_col_part_ = 0;
  if(ref_row_part_.has_ref_set())		ref_row_part_.obj().resize(0);
  if(ref_col_part_.has_ref_set())		ref_col_part_.obj().resize(0);
  if(ref_ele_.has_ref_set())			ref_ele_.obj().resize(0,element_type());
  if(ref_part_start_.has_ref_set())	ref_part_start_.obj().resize(0);
}

template <class T_Indice, class T_Value>
void COOMatrixPartitionedView<T_Indice,T_Value>::get_row_part(indice_type row_part[]) const
{
  assert_initialized();
  const std::vector<size_type> &_row_part = ref_row_part_.const_obj();
  std::copy(_row_part.begin(), _row_part.end(), row_part);
}

template <class T_Indice, class T_Value>
void COOMatrixPartitionedView<T_Indice,T_Value>::get_col_part(indice_type col_part[]) const
{
  assert_initialized();
  const std::vector<size_type> &_col_part = ref_col_part_.const_obj();
  std::copy(_col_part.begin(), _col_part.end(), col_part);
}

template <class T_Indice, class T_Value>
COOMatrixPartitionedView<T_Indice,T_Value>::size_type
COOMatrixPartitionedView<T_Indice,T_Value>::part_num(const vector_size_type& part
  , size_type indice)
{
  return std::upper_bound(part.begin(),part.end(),indice) - part.begin();
}

template <class T_Indice, class T_Value>
COOMatrixPartitionedView<T_Indice,T_Value>::partition_type
COOMatrixPartitionedView<T_Indice,T_Value>::create_partition(Range1D rng_overall_p) const
{
  assert_initialized();
  // get reference to data structures
  const std::vector<size_type>
    &row_part	= ref_row_part_.const_obj(),
    &col_part	= ref_col_part_.const_obj(),
    &part_start	= ref_part_start_.const_obj();
  ele_type
    &_ele		= const_cast<ref_ele_type&>(ref_ele_).obj(); // This is ok
  // Get upper and lower overall, row and column partition numbers
  rng_overall_p = DenseLinAlgPack::full_range(rng_overall_p,1,num_row_part_*num_col_part_);
  size_type	l_p		= rng_overall_p.lbound(),
        u_p		= rng_overall_p.ubound(),
        l_r_p	= imp_row_part_num(l_p),
        l_c_p	= imp_col_part_num(l_p),
        u_r_p	= imp_row_part_num(u_p),
        u_c_p	= imp_col_part_num(u_p);
  // build argument list for creation of the partition.
  size_type
    rows	= row_part[u_r_p] - row_part[l_r_p - 1],
    cols	= col_part[u_c_p] - col_part[l_c_p - 1],
    nz		= part_start[u_p] - part_start[l_p - 1];
  element_type
    *ele		= &_ele[0] + part_start[l_p - 1];
  difference_type
    row_offset	= - (row_part[l_r_p - 1] - 1),
    col_offset	= - (col_part[l_c_p - 1] - 1);

  return partition_type(rows,cols,nz,ele,row_offset,col_offset);
}

} // end namespace AbstractLinAlgPack 

#endif // COO_MATRIX_PARTITIONED_VIEW_CLASS_DEF_H
