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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_InterlacedOperator.hpp"

#include "Epetra_Map.h"

#include "EpetraExt_BlockUtility.h"

Stokhos::InterlacedOperator::
InterlacedOperator(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_CrsGraph>& determ_graph,
  const Teuchos::RCP<Teuchos::ParameterList>& params) : 
  EpetraExt::BlockCrsMatrix(*(epetraCijk_->getStochasticGraph()), 
                            *determ_graph,
			    *sg_comm_),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  Cijk(epetraCijk->getParallelCijk()),
  block_ops(),
  scale_op(true),
  include_mean(true),
  only_use_linear(false),
  determOffset_(EpetraExt::BlockUtility::CalculateOffset(determ_graph->RowMap()))
{
  scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
  include_mean = params->get("Include Mean", true);
  only_use_linear = params->get("Only Use Linear Terms", false);
}

Stokhos::InterlacedOperator::
~InterlacedOperator()
{
}

void 
Stokhos::InterlacedOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops)
{
  block_ops = ops;

  // Zero out matrix
  this->PutScalar(0.0);

  // Compute loop bounds
  Cijk_type::k_iterator k_begin = Cijk->k_begin();
  Cijk_type::k_iterator k_end = Cijk->k_end();
  if (!include_mean && index(k_begin) == 0)
    ++k_begin;
  if (only_use_linear) {
    int dim = sg_basis->dimension();
    k_end = Cijk->find_k(dim+1);
  }

  // Assemble matrix
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int k = index(k_it);
    Teuchos::RCP<Epetra_RowMatrix> block = 
      Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(block_ops->getCoeffPtr(k), 
						  true);
    for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	 j_it != Cijk->j_end(k_it); ++j_it) {
      int j = epetraCijk->GCID(index(j_it));
      for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	   i_it != Cijk->i_end(j_it); ++i_it) {
	int i = epetraCijk->GRID(index(i_it));
	double c = value(i_it);
	if (scale_op)
	  c /= norms[i];
	this->SumIntoGlobalBlock_Deterministic(c, *block, i, j);
      }
    }
  }
}

Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > 
Stokhos::InterlacedOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly > 
Stokhos::InterlacedOperator::
getSGPolynomial() const
{
  return block_ops;
}

void Stokhos::InterlacedOperator::
SumIntoGlobalBlock_Deterministic(double alpha,const Epetra_RowMatrix & determBlock,int Row,int Col)
{
  const Epetra_Map & determMap = determBlock.RowMatrixRowMap();
  const Epetra_Map & determColMap = determBlock.RowMatrixColMap();

  // This routine copies entries of a BaseMatrix into big  BlockCrsMatrix
  // It performs the following operation on the global IDs row-by-row
  // this->val[i+rowOffset][j+ColOffset] = BaseMatrix.val[i][j]

  int MaxIndices = determBlock.MaxNumEntries();
  std::vector<int> Indices(MaxIndices);
  std::vector<double> Values(MaxIndices);
  int NumIndices;
  int ierr=0;

  for (int i=0; i<determMap.NumMyElements(); i++) {
    determBlock.ExtractMyRowCopy( i, MaxIndices, NumIndices, &Values[0], &Indices[0] );

    // Convert to BlockMatrix Global numbering scheme
    for( int l = 0; l < NumIndices; ++l ) {
       Indices[l] = Col +  COffset_*determColMap.GID(Indices[l]);
       Values[l] *= alpha;
    }

    int BaseRow = determMap.GID(i);
    ierr = this->SumIntoGlobalValues(ROffset_*BaseRow + Row, NumIndices, &Values[0], &Indices[0]);
    if (ierr != 0) cout << "WARNING InterlacedOperator::SumIntoBlock_Deterministic SumIntoGlobalValues err = " << ierr <<
            "\n\t  Row " << ROffset_*BaseRow + Row << ", Col start " << Indices[0] << endl;

  }
}
