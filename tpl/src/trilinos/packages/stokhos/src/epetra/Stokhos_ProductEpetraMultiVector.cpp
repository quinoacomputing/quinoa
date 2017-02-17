// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector() :
  ProductContainer<Epetra_MultiVector>()
{
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(const Teuchos::RCP<const Epetra_BlockMap>& block_map) :
  ProductContainer<Epetra_MultiVector>(block_map)
{
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
  int num_vectors) : 
  ProductContainer<Epetra_MultiVector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*coeff_map,*block_map, *product_comm))),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map,
						  num_vectors)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
  int num_vectors) :
  ProductContainer<Epetra_MultiVector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(product_map_),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map,
						  num_vectors)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(
  const Teuchos::RCP<const Epetra_BlockMap>& block_map,
  const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
  const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
  Epetra_DataAccess CV,
  const Epetra_MultiVector& block_vector) :
  ProductContainer<Epetra_MultiVector>(block_map),
  coeff_map(coeff_map_),
  product_comm(product_comm_),
  product_map(product_map_),
  bv(Teuchos::rcp(new EpetraExt::BlockMultiVector(CV, *coeff_map, 
						  block_vector)))
{
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}
    
Stokhos::ProductEpetraMultiVector::
ProductEpetraMultiVector(const Stokhos::ProductEpetraMultiVector& v) :
  ProductContainer<Epetra_MultiVector>(v),
  coeff_map(v.coeff_map),
  product_comm(v.product_comm),
  product_map(v.product_map),
  bv(v.bv)
{
}

Stokhos::ProductEpetraMultiVector::
~ProductEpetraMultiVector() {}

Stokhos::ProductEpetraMultiVector& 
Stokhos::ProductEpetraMultiVector::
operator=(const Stokhos::ProductEpetraMultiVector& v) {
  ProductContainer<Epetra_MultiVector>::operator=(v);
  coeff_map = v.coeff_map;
  product_comm = v.product_comm;
  product_map = v.product_map;
  bv = v.bv;  // Note this is a shallow copy, which is consistent with above
  return *this;
}

Stokhos::ProductEpetraMultiVector& 
Stokhos::ProductEpetraMultiVector::
operator=(const Epetra_MultiVector& v) {
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
  return *this;
}

void 
Stokhos::ProductEpetraMultiVector::
assignToBlockMultiVector(Epetra_MultiVector& v) const 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      v.Update(1.0, *bv, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(block_v.GetBlock(i)) = *(coeff_[i]);
    }
  }
}

void 
Stokhos::ProductEpetraMultiVector::
assignFromBlockMultiVector(const Epetra_MultiVector& v) 
{
  if (this->size() > 0) {
    if (bv != Teuchos::null)
      bv->Update(1.0, v, 0.0);
    else {
      EpetraExt::BlockMultiVector block_v(View, *coeff_map, v);
      for (int i=0; i<this->size(); i++)
	*(coeff_[i]) = *(block_v.GetBlock(i));
    }
  }
}

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::ProductEpetraMultiVector::
coefficientMap() const {
  return coeff_map;
}

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::ProductEpetraMultiVector::
productMap() const {
  return product_map;
}

Teuchos::RCP<const EpetraExt::MultiComm> 
Stokhos::ProductEpetraMultiVector::
productComm() const {
  return product_comm;
}

int
Stokhos::ProductEpetraMultiVector::
numVectors() const {
  if (bv != Teuchos::null)
    return bv->NumVectors();
  else if (this->size() > 0 && this->coeff_[0] != Teuchos::null)
    return this->coeff_[0]->NumVectors();
  return -1;
}
      
void 
Stokhos::ProductEpetraMultiVector::
reset(const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
      int num_vectors) 
{
  ProductContainer<Epetra_MultiVector>::reset(block_map);
  coeff_map = coeff_map_;
  product_comm = product_comm_;
  product_map = 
    Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*coeff_map,
							   *block_map, 
							   *product_comm));
  bv = Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map, 
						    num_vectors));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

void 
Stokhos::ProductEpetraMultiVector::
reset(const Teuchos::RCP<const Epetra_BlockMap>& block_map,
      const Teuchos::RCP<const Epetra_BlockMap>& coeff_map_,
      const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
      const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
      int num_vectors) 
{
  ProductContainer<Epetra_MultiVector>::reset(block_map);
  coeff_map = coeff_map_;
  product_comm = product_comm_;
  product_map = product_map_;
  bv = Teuchos::rcp(new EpetraExt::BlockMultiVector(*coeff_map, *product_map,
						    num_vectors));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

void 
Stokhos::ProductEpetraMultiVector::
resetCoefficients(Epetra_DataAccess CV,
		  const Epetra_MultiVector& block_vector) 
{
  bv = 
    Teuchos::rcp(new EpetraExt::BlockMultiVector(CV, *coeff_map, block_vector));
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}

Teuchos::RCP<EpetraExt::BlockMultiVector> 
Stokhos::ProductEpetraMultiVector::
getBlockMultiVector() 
{
  return bv;
}

Teuchos::RCP<const EpetraExt::BlockMultiVector> 
Stokhos::ProductEpetraMultiVector::
getBlockMultiVector() const 
{
  return bv;
}

void 
Stokhos::ProductEpetraMultiVector::
setBlockMultiVector(const Teuchos::RCP<EpetraExt::BlockMultiVector>& block_vec) 
{
  bv = block_vec;
  for (int i=0; i<this->size(); i++)
    this->setCoeffPtr(i, bv->GetBlock(i));
}
