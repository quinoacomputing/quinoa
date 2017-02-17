// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SPARSE_3_TENSOR_UTILITIES_HPP
#define STOKHOS_SPARSE_3_TENSOR_UTILITIES_HPP

#include "Stokhos_Sparse3Tensor.hpp"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_BlockMap.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

namespace Stokhos {

  //! Build an Epetra_CrsGraph from a sparse 3 tensor
  /*!
   * Builds a sparse graph from a sparse 3 tensor by summing over the third
   * index.  This graph then represents the sparsity pattern of the stochastic
   * part of the block stochastic Galerkin operator.  Redistributing the graph
   * should then provide a suitable parallel distribution for block
   * stochastic Galerkin linear solves.
   */
  template <typename ordinal_type, typename value_type>
  Teuchos::RCP<Epetra_CrsGraph>
  sparse3Tensor2CrsGraph(
    const Stokhos::OrthogPolyBasis<ordinal_type,value_type>& basis,
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_Comm& comm) 
  {
    typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk_type;

    // Number of stochastic rows
    ordinal_type num_rows = basis.size();

    // Replicated local map
    Epetra_LocalMap map(num_rows, 0, comm);

    // Graph to be created
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Teuchos::rcp(new Epetra_CrsGraph(Copy, map, 0));
    
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    for (typename Cijk_type::k_iterator k_it=Cijk.k_begin(); 
	 k_it!=Cijk.k_end(); ++k_it) {
      for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	   j_it != Cijk.j_end(k_it); ++j_it) {
	ordinal_type j = index(j_it);
	for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it);
	     i_it != Cijk.i_end(j_it); ++i_it) {
	  ordinal_type i = index(i_it);
	  graph->InsertGlobalIndices(i, 1, &j);
	}
      }
    }

    // Sort, remove redundencies, transform to local, ...
    graph->FillComplete();

    return graph;
  }

  //! Build an Epetra_CrsGraph from a sparse 3 tensor
  /*!
   * Builds a sparse graph from a sparse 3 tensor by summing over the third
   * index.  This graph then represents the sparsity pattern of the stochastic
   * part of the block stochastic Galerkin operator.  Redistributing the graph
   * should then provide a suitable parallel distribution for block
   * stochastic Galerkin linear solves.
   */
  template <typename ordinal_type, typename value_type>
  Teuchos::RCP<Epetra_CrsGraph>
  sparse3Tensor2CrsGraph(
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_BlockMap& map) 
  {
    typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk_type;

    // Graph to be created
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Teuchos::rcp(new Epetra_CrsGraph(Copy, map, 0));
    
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    for (typename Cijk_type::k_iterator k_it=Cijk.k_begin(); 
	 k_it!=Cijk.k_end(); ++k_it) {
      for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	   j_it != Cijk.j_end(k_it); ++j_it) {
	ordinal_type j = index(j_it);
	for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it);
	     i_it != Cijk.i_end(j_it); ++i_it) {
	  ordinal_type i = index(i_it);
	  graph->InsertGlobalIndices(i, 1, &j);
	}
      }
    }

    // Sort, remove redundencies, transform to local, ...
    graph->FillComplete();

    return graph;
  }

  template <typename ordinal_type, typename value_type>
  void
  sparse3Tensor2MatrixMarket(
    const Stokhos::OrthogPolyBasis<ordinal_type,value_type>& basis,
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_Comm& comm,
    const std::string& file)
  {
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Stokhos::sparse3Tensor2CrsGraph(basis, Cijk, comm);
    Epetra_CrsMatrix mat(Copy, *graph);
    mat.FillComplete();
    mat.PutScalar(1.0);
    EpetraExt::RowMatrixToMatrixMarketFile(file.c_str(), mat);
  }

  template <typename ordinal_type, typename value_type>
  void
  sparse3Tensor2MatrixMarket(
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_BlockMap& map,
    const std::string& file)
  {
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Stokhos::sparse3Tensor2CrsGraph(Cijk, map);
    Epetra_CrsMatrix mat(Copy, *graph);
    mat.FillComplete();
    mat.PutScalar(1.0);
    EpetraExt::RowMatrixToMatrixMarketFile(file.c_str(), mat);
  }

} // namespace Stokhos

#endif // SPARSE_3_TENSOR_UTILITIES_HPP
