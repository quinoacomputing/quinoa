/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_STOKHOSBASISWRAPPER_H
#define SUNDANCE_STOKHOSBASISWRAPPER_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMap.hpp"
#include "SundanceSpectralBasisBase.hpp"

#ifdef HAVE_SUNDANCE_STOKHOS

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"

#include "cijk.h"
#include "chaos.h"



namespace Sundance
{
using Teuchos::RCP;

/** 
 * 
 */
class  StokhosBasisWrapper : public SpectralBasisBase
{
public: 
  /** */
  typedef Stokhos::OrthogPolyBasis<int, double> PolyBasis;
  /** */
  typedef Stokhos::OneDOrthogPolyBasis<int, double> PolyBasis1D;

  /** Construct a basis */
  StokhosBasisWrapper(const RCP<const PolyBasis>& basis);
  /** Construct a basis */
  StokhosBasisWrapper(const RCP<const PolyBasis1D>& basis);
    
  /** Return the dim of the Spectral Basis */
  int getDim() const {return basis_->dimension();}
  
  /** Return the order of the Spectral Basis */
  int getOrder() const {return basis_->order();}
  
  /** Return the maximum number of terms */
  int nterms() const {return basis_->size();}
  
  /** Return the basis element stored in the basis array index */
  int getElement(int i) const {return i;}
  
  /** expectation operator */
  double expectation(int i, int j, int k);
  
  /** Write to a std::string */
  std::string toString() const {return basis_->getName();}
  
  /* */
  GET_RCP(SpectralBasisBase);
  
  /** Ordering operator */
  virtual bool lessThan(const SpectralBasisBase* other) const ;

private:
  RCP<const PolyBasis> basis_;
  RCP<Stokhos::Sparse3Tensor<int, double> > cijk_;
  
  /** Fill the c_ijk array */
  void fillCijk() ;
};
}

#endif
#endif
