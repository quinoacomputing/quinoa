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

#include "Stokhos_BlockDiagonalOperator.hpp"
#include "EpetraExt_BlockMultiVector.h"
#include "Epetra_Map.h"

Stokhos::BlockDiagonalOperator::
BlockDiagonalOperator(
  const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm_,
  int num_mp_blocks_,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
  const Teuchos::RCP<const Epetra_Map>& range_base_map_,
  const Teuchos::RCP<const Epetra_Map>& domain_mp_map_,
  const Teuchos::RCP<const Epetra_Map>& range_mp_map_) : 
  label("Stokhos Block Diagonal Operator"),
  mp_comm(mp_comm_),
  num_mp_blocks(num_mp_blocks_),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  domain_mp_map(domain_mp_map_),
  range_mp_map(range_mp_map_),
  block_ops(),
  useTranspose(false)
{
}

Stokhos::BlockDiagonalOperator::
~BlockDiagonalOperator()
{
}

void 
Stokhos::BlockDiagonalOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::ProductEpetraOperator >& ops)
{
  block_ops = ops;
}

Teuchos::RCP< Stokhos::ProductEpetraOperator > 
Stokhos::BlockDiagonalOperator::
getMPOps()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::ProductEpetraOperator > 
Stokhos::BlockDiagonalOperator::
getMPOps() const
{
  return block_ops;
}

int 
Stokhos::BlockDiagonalOperator::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (int i=0; i<num_mp_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::BlockDiagonalOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  Teuchos::RCP<const Epetra_BlockMap> input_base_map, result_base_map;
  if (useTranspose) {
    input_base_map = range_base_map;
    result_base_map = domain_base_map;
  }
  else {
    input_base_map = domain_base_map;
    result_base_map = range_base_map;
  }
  EpetraExt::BlockMultiVector mp_input(View, *input_base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *result_base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    (*block_ops)[i].Apply(*(mp_input.GetBlock(i)), *(mp_result.GetBlock(i)));
  }

  return 0;
}

int 
Stokhos::BlockDiagonalOperator::ApplyInverse(const Epetra_MultiVector& Input, 
					     Epetra_MultiVector& Result) const
{
  Teuchos::RCP<const Epetra_BlockMap> input_base_map, result_base_map;
  if (useTranspose) {
    input_base_map = domain_base_map;
    result_base_map = range_base_map;
  }
  else {
    input_base_map = range_base_map;
    result_base_map = domain_base_map;
  }
  EpetraExt::BlockMultiVector mp_input(View, *input_base_map, Input);
  EpetraExt::BlockMultiVector mp_result(View, *range_base_map, Result);
  for (int i=0; i<num_mp_blocks; i++) {
    (*block_ops)[i].ApplyInverse(*(mp_input.GetBlock(i)), 
				 *(mp_result.GetBlock(i)));
  }

  return 0;
}

double 
Stokhos::BlockDiagonalOperator::NormInf() const
{
  double product_nrm = 0.0;
  for (int i=0; i<num_mp_blocks; i++) {
    double nrm = (*block_ops)[i].NormInf();
    if (nrm > product_nrm)
      product_nrm = nrm;
  }

  return product_nrm;
}


const char* 
Stokhos::BlockDiagonalOperator::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::BlockDiagonalOperator::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::BlockDiagonalOperator::HasNormInf() const
{
  if (num_mp_blocks == 0)
    return false;
  return (*block_ops)[0].HasNormInf();
}

const Epetra_Comm & 
Stokhos::BlockDiagonalOperator::Comm() const
{
  return *mp_comm;
}
const Epetra_Map& 
Stokhos::BlockDiagonalOperator::OperatorDomainMap() const
{
  if (useTranspose)
    return *range_mp_map;
  return *domain_mp_map;
}

const Epetra_Map& 
Stokhos::BlockDiagonalOperator::OperatorRangeMap() const
{
  if (useTranspose)
    return *domain_mp_map;
  return *range_mp_map;
}
