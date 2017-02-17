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

#include "Stokhos_AdaptivityManager.hpp"
#include "Stokhos_AdaptivityUtils.hpp"
#include "Stokhos_BasisInteractionGraph.hpp"

#include "EpetraExt_BlockVector.h"
#include "EpetraExt_RowMatrixOut.h"

#ifdef HAVE_STOKHOS_BOOST
#include <boost/functional/hash.hpp>
#endif

#ifdef HAVE_STOKHOS_BOOST // we have boost, use the hash Stokhos, use the hash!
std::size_t Stokhos::AdaptivityManager::Sparse3TensorHash::IJKHash::
operator()(const Stokhos::AdaptivityManager::Sparse3TensorHash::IJK & ijk) const
{
   std::size_t seed=0;
   boost::hash_combine(seed,ijk.i_);
   boost::hash_combine(seed,ijk.j_);
   boost::hash_combine(seed,ijk.k_);
   return seed;
}

Stokhos::AdaptivityManager::Sparse3TensorHash::Sparse3TensorHash(const Stokhos::Sparse3Tensor<int,double> & Cijk)
{
   typedef Stokhos::Sparse3Tensor<int,double>::k_iterator k_iterator;
   typedef Stokhos::Sparse3Tensor<int,double>::kj_iterator kj_iterator;
   typedef Stokhos::Sparse3Tensor<int,double>::kji_iterator kji_iterator;

   for(k_iterator k_it = Cijk.k_begin();k_it!=Cijk.k_end();k_it++) {
      int k = *k_it;
      for(kj_iterator j_it = Cijk.j_begin(k_it);j_it!=Cijk.j_end(k_it);j_it++) {
         int j = *j_it;
         for(kji_iterator i_it = Cijk.i_begin(j_it);i_it!=Cijk.i_end(j_it);i_it++) {
            int i = *i_it;
            hashMap_[IJK(i,j,k)] = i_it.value();
         }
      }
   }
}

double Stokhos::AdaptivityManager::Sparse3TensorHash::getValue(int i,int j,int k) const
{
   boost::unordered_map<IJK,double>::const_iterator itr;
   itr = hashMap_.find(IJK(i,j,k));

   if(itr==hashMap_.end()) return 0.0;

   return itr->second;
}
#else // no BOOST, just default to the slow thing
Stokhos::AdaptivityManager::Sparse3TensorHash::Sparse3TensorHash(const Stokhos::Sparse3Tensor<int,double> & Cijk)
      : Cijk_(Cijk)
{ }

double Stokhos::AdaptivityManager::Sparse3TensorHash::getValue(int i,int j,int k) const
{
   return Cijk_.getValue(i,j,k);
}
#endif

Stokhos::AdaptivityManager::AdaptivityManager(
         const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& sg_master_basis,
         const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & sg_basis_row_dof,
         const Epetra_CrsGraph & determ_graph,
         bool onlyUseLinear,int kExpOrder,
         bool scaleOp)
   : sg_master_basis_(sg_master_basis), sg_basis_row_dof_(sg_basis_row_dof), scaleOp_(scaleOp)
{
    rowMap_ = adapt_utils::buildAdaptedRowMapAndOffsets(determ_graph.Comm(),sg_basis_row_dof_,myRowGidOffsets_);

    setupWithGraph(determ_graph,onlyUseLinear,kExpOrder);
}

Stokhos::AdaptivityManager::AdaptivityManager(
         const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& sg_master_basis,
         const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & sg_basis_row_dof,
         const Epetra_Comm & comm,
         bool scaleOp)
   : sg_master_basis_(sg_master_basis), sg_basis_row_dof_(sg_basis_row_dof), scaleOp_(scaleOp)
{
   rowMap_ = adapt_utils::buildAdaptedRowMapAndOffsets(comm,sg_basis_row_dof_,myRowGidOffsets_);
}

Teuchos::RCP<Epetra_CrsMatrix> 
Stokhos::AdaptivityManager::buildMatrixFromGraph() const
{
   return Teuchos::rcp(new Epetra_CrsMatrix(Copy,*graph_));
}

void Stokhos::AdaptivityManager::setupWithGraph(const Epetra_CrsGraph & determGraph,bool onlyUseLinear,int kExpOrder) 
{
   graph_ = adapt_utils::buildAdaptedGraph(determGraph, sg_master_basis_, sg_basis_row_dof_, onlyUseLinear, kExpOrder);

   adapt_utils::buildAdaptedColOffsets(determGraph,myRowGidOffsets_,myColGidOffsets_);
   adapt_utils::buildColBasisFunctions(determGraph,sg_master_basis_,sg_basis_row_dof_,sg_basis_col_dof_);
}

/** Setup operator
  */
void 
Stokhos::AdaptivityManager::
setupOperator(Epetra_CrsMatrix & A,const Sparse3Tensor<int,double> & Cijk,Stokhos::EpetraOperatorOrthogPoly & poly,
              bool onlyUseLinear,bool includeMean) const
{
   typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

   // build the sparse hash only once
   Sparse3TensorHash hashLookup(Cijk);

   // Zero out matrix
   A.PutScalar(0.0);

   // Compute loop bounds
   Cijk_type::k_iterator k_begin = Cijk.k_begin();
   Cijk_type::k_iterator k_end = Cijk.k_end();
   if (!includeMean && index(k_begin) == 0)
     ++k_begin;
   if (onlyUseLinear) {
     int dim = sg_master_basis_->dimension();
     k_end = Cijk.find_k(dim+1);
   }
 
   // Assemble matrix
   // const Teuchos::Array<double>& norms = sg_basis->norm_squared();
   for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int k = index(k_it);
      Teuchos::RCP<Epetra_CrsMatrix> block = 
        Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(poly.getCoeffPtr(k), 
 						    true);

      // add in matrix k
      sumInOperator(A,hashLookup,k,*block);
   }
}

void
Stokhos::AdaptivityManager::
sumInOperator(Epetra_CrsMatrix & A,const Stokhos::Sparse3Tensor<int,double> & Cijk,int k,const Epetra_CrsMatrix & J_k) const
{
   // this allows the simple interface of taking a Sparse3Tensor but immediately computes
   // the sparse hash (if boost is enabled)
  
   Sparse3TensorHash hashLookup(Cijk);
   sumInOperator(A,hashLookup,k,J_k);
}

void
Stokhos::AdaptivityManager::
sumInOperator(Epetra_CrsMatrix & A,const Stokhos::AdaptivityManager::Sparse3TensorHash & Cijk,int k,const Epetra_CrsMatrix & J_k) const
{
   TEUCHOS_ASSERT(J_k.NumMyRows() == int(sg_basis_row_dof_.size()));
   TEUCHOS_ASSERT(J_k.NumMyCols() == int(sg_basis_col_dof_.size()));

   const Teuchos::Array<double> & normValues = sg_master_basis_->norm_squared();

   // loop over deterministic rows 
   for(int localM=0;localM<J_k.NumMyRows();localM++) {
     // int m = J_k.GRID(localM); // unused

      // grab row basis
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > rowStochBasis 
            = sg_basis_row_dof_[localM]; 
 
      // grab row from deterministic system
      int d_numEntries;
      int * d_Indices;
      double * d_Values;
     
      J_k.ExtractMyRowView(localM,d_numEntries,d_Values,d_Indices);
      
      // loop over stochastic degrees of freedom of this row
      for(int rb_i=0;rb_i<rowStochBasis->size();rb_i++) {
         int i = sg_master_basis_->index(rowStochBasis->term(rb_i));

         double normValue = normValues[i]; // sg_master_basis->norm_squared(i);
         
         int sg_m = getGlobalRowId(localM,rb_i);

         // we wipe out old values, capacity should gurantee
         // we don't allocate more often than neccessary!
         std::vector<int> sg_indices;
         std::vector<double> sg_values;

         // sg_indices.resize(0); 
         // sg_values.resize(0);

         // loop over each column
         for(int colInd=0;colInd<d_numEntries;colInd++) {
            int localN = d_Indices[colInd];  // grab local deterministic column id

            // grab row basis
            Teuchos::RCP<const Stokhos::ProductBasis<int,double> > colStochBasis 
                  = sg_basis_col_dof_[localN]; 

            // build values array
            for(int cb_j=0;cb_j<colStochBasis->size();cb_j++) {
               int j = sg_master_basis_->index(colStochBasis->term(cb_j));
               int sg_n = getGlobalColId(localN,cb_j);
               double cijk = Cijk.getValue(i,j,k); 

               // no reason to work it in!
               if(cijk==0) continue;

               if(scaleOp_)
                  cijk = cijk/normValue;

               sg_indices.push_back(sg_n);
               sg_values.push_back(cijk*d_Values[colInd]);
            }
         }

         // add in matrix values
         A.SumIntoGlobalValues(sg_m,sg_indices.size(),&sg_values[0],&sg_indices[0]);
      }
   }
}

/** Copy to an adaptive vector from a set of blocked vectors
  */
void Stokhos::AdaptivityManager::copyToAdaptiveVector(const Stokhos::EpetraVectorOrthogPoly & x_sg,Epetra_Vector & x) const
{
   Teuchos::RCP<const EpetraExt::BlockVector> x_sg_bv = x_sg.getBlockVector();

   // copy from adapted vector to deterministic
   for(std::size_t i=0;i<sg_basis_row_dof_.size();i++) {
      int P_i = getRowStochasticBasisSize(i); 
      int localId = rowMap_->LID(getGlobalRowId(i,0));

      for(int j=0;j<P_i;j++,localId++) {
         int blk = sg_master_basis_->index(sg_basis_row_dof_[i]->term(j));
         x[localId] = x_sg_bv->GetBlock(blk)->operator[](i);
      }
   }
}

/** Copy from an adaptive vector to a set of blocked vectors
  */
void Stokhos::AdaptivityManager::copyFromAdaptiveVector(const Epetra_Vector & x,Stokhos::EpetraVectorOrthogPoly & x_sg) const
{
   int numBlocks = x_sg.size();
   Teuchos::RCP<EpetraExt::BlockVector> x_sg_bv = x_sg.getBlockVector();

   // zero out determinstic vectors
   for(int blk=0;blk<numBlocks;blk++)
      x_sg_bv->GetBlock(blk)->PutScalar(0.0);

   // copy from adapted vector to deterministic
   for(std::size_t i=0;i<sg_basis_row_dof_.size();i++) {
      int P_i = getRowStochasticBasisSize(i); 
      int localId = rowMap_->LID(getGlobalRowId(i,0));

      for(int j=0;j<P_i;j++,localId++) {
         int blk = sg_master_basis_->index(sg_basis_row_dof_[i]->term(j));
         x_sg_bv->GetBlock(blk)->operator[](i) = x[localId];
      }
   }
}
