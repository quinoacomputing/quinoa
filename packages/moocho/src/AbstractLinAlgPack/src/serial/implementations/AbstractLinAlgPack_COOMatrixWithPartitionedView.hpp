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

#ifndef COO_MATRIX_WITH_PARTITIONED_VIEW_H
#define COO_MATRIX_WITH_PARTITIONED_VIEW_H

#include "AbstractLinAlgPack_COOMatrixClass.hpp"
#include "AbstractLinAlgPack_COOMatrixPartitionedViewClass.hpp"

namespace AbstractLinAlgPack {

/** \brief Aggregation of a COO matrix and a partitioned view of it.
  *
  * This class represents the aggregation of a COOMatrix and a COOMatrixParitionedView.
  * This class is designed to help avoid mistakes that may happen when the underlying
  * COO matrix is modified and the partitioned view becomes obsolete.  Therefore, operations
  * that may make the partitioned view obsolete are encapsulated to delete the partitioned view.
  *
  * Therefore, non-const references to the underlying COOMatrix and COOMatrixPartitioned view
  * are not provided.  Analogs for the non-const member functions
  * are provided and keep track of the book keeping for you.  This is a very light weight class
  * in terms of overhead.  The only exception to this are the rows() and cols() member functions.
  * They are included to allow for use with the MatrixWithOpConcreteEncap<M> class.
  *
  * The const interfaces to these objects can be accessed using the coom() and coom_view()
  * member functions.
  */
class COOMatrixWithPartitionedView {
public:

  // /////////////////////////////////////////////////////////////
  /** @name Public Types */
  //@{

  /** \brief . */
  typedef COOMatrixPartitionedView<indice_type,value_type>	partitioned_view_type;
  /** \brief . */
  typedef partitioned_view_type::partition_type				partition_type;
  /** \brief . */
  typedef partitioned_view_type::EPartitionOrder				EPartitionOrder;

  //@}

  /// Calls coom_view().rows() if the partitioned view has been initialize and coom().rows() if not
  size_type rows() const {
    return coom_view_.is_initialized() ? coom_view_.rows() : coom_.rows();
  }

  /// Calls coom_view().cols() if the partitioned view has been initialize and coom().cols() if not
  size_type cols() const {
    return coom_view_.is_initialized() ? coom_view_.cols() : coom_.cols();
  }

  /// Assignment operator
  COOMatrixWithPartitionedView& operator=(const COOMatrixWithPartitionedView& m) {
    coom_ = m.coom_;
    coom_view_.bind(coom_view_);
    return *this;
  }
  
  /** \brief Allow assignment to a COOMatrix.
    *
    * After the assignment #this# will not have its partitioned view initialized.
    */
  COOMatrixWithPartitionedView& operator=(const COOMatrix& m) {
    coom_ = m;
    coom_view_.free();
    return *this;
  }
  
  // /////////////////////////////////////////////////////////////
  /** @name COOMatrix non-const encapsulated interface.
    *
    * The default constructor and copy constructor are allowed.
    */
  //@{

  /** \brief Resize for a #rows# by #cols# sparse matrix with #nz# elements.
    *
    * If there was a partitioned view set up then this will destory it.
    */
  void resize(size_type rows, size_type cols, size_type nz) {
    coom_view_.free();
    coom_.resize(rows,cols,nz);
  }

  /** \brief Return pointer to raw storage array (length #nz()#) for the values of the non-zero elements.
    *
    * This operation will not result in a loss of the partitioned view.
    */
  value_type*							val() {
    return coom_.val();
  }
  /** \brief Return pointer to raw storage array (length #nz()#) for the row indices of the non-zero elements
    *
    * This operation will result in a loss of the partitioned view.
    */
  indice_type*						ivect() {
    coom_view_.free();
    return coom_.ivect();
  }
  /** \brief Return pointer to raw storage array (length #nz()#) for the column indices of the non-zero elements
    *
    * This operation will result in a loss of the partitioned view.
    */
  indice_type*						jvect() {
    coom_view_.free();
    return coom_.jvect();
  }

  /** \brief Initialize from an input stream.
    *
    * This operation will result in a loss of the partitioned view.
    */
  void initialize(std::istream& istrm) {
    coom_view_.free();
    coom_.initialize(istrm);
  }

  //@}

  /// Return a const referece to the COOMatrix
  const COOMatrix& coom() const {
    return coom_;
  }

  // ////////////////////////////////////////////////////////////////
  /** @name Non-const COOMatrixPartitionedView interface */
  //@{

  /** \brief Crete a view to the COO matrix.
    *
    * Calls create_view on the partitioned view object using the data
    * from the COOMatrix object.
    */
  void create_view(
        const size_type		row_perm[]
      , const size_type		col_perm[]
      , const size_type		num_row_part
      , const size_type		row_part[]
      , const size_type		num_col_part
      , const size_type		col_part[]
      , const EPartitionOrder	partition_order )
  {
    coom_view_.create_view(coom_.rows(),coom_.cols(),coom_.nz(),coom_.val(),coom_.const_ivect()
      ,coom_.const_jvect(),row_perm,col_perm,num_row_part,row_part,num_col_part,col_part
      ,partition_order);
  }

  /** \brief . */
  partition_type partition(size_type overall_p) {
    return coom_view_.partition(overall_p);
  }

  /** \brief . */
  partition_type partition(size_type row_p, size_type col_p) {
    return coom_view_.partition(row_p, col_p);
  }

  /** \brief . */
  partition_type partition(Range1D rng_overall_p) {
    return coom_view_.partition(rng_overall_p);
  }

  //@}

  /// Return a const referece to the COOMatrixPartitionedView object
  const partitioned_view_type& coom_view() const {
    return coom_view_;
  }


private:
  COOMatrix				coom_;
  partitioned_view_type	coom_view_;


};	// end class COOMatrixWithPartitionedView


} // end namespace AbstractLinAlgPack

#endif // COO_MATRIX_WITH_PARTITIONED_VIEW_H
