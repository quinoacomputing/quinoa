// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_EXAMPLE_MYCELL_HPP
#define PHX_EXAMPLE_MYCELL_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Dimension.hpp"
#include "Shards_Array.hpp"

class MyCell {
  
public:

  MyCell();
  
  virtual ~MyCell() {}
  
  shards::Array<double,shards::NaturalOrder,Node,Dim>& getNodeCoordinates();
  
  shards::Array<double,shards::NaturalOrder,QuadPoint,Node>& 
  getBasisFunctions();
  
  shards::Array<double,shards::NaturalOrder,QuadPoint,Node,Dim>& 
  getBasisFunctionGradients();
  
  std::size_t localIndex();

  void setLocalIndex(std::size_t index);

private:
  
  std::size_t local_index_;

  Teuchos::ArrayRCP<double> m_coords_mem;
  
  Teuchos::ArrayRCP<double> m_phi_mem;
  
  Teuchos::ArrayRCP<double> m_grad_phi_mem;

  shards::Array<double,shards::NaturalOrder,Node,Dim> m_coords;
  
  shards::Array<double,shards::NaturalOrder,QuadPoint,Node> m_phi;

  shards::Array<double,shards::NaturalOrder,QuadPoint,Node,Dim> m_grad_phi;

};

#endif
