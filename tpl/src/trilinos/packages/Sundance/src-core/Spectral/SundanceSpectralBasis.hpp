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

#ifndef SUNDANCE_SPECTRALBASIS_H
#define SUNDANCE_SPECTRALBASIS_H

#include "SundanceDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaHandleable.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceSpectralBasisBase.hpp"
#include "SundanceStokhosBasisWrapper.hpp"



namespace Sundance
{
using Playa::Handle;
using Playa::Handleable;



/** Doxygen doc for SpectralBasis */
class SpectralBasis : public Playa::Handle<SpectralBasisBase>
{
public:

  /* boilerplate handle ctors */
  HANDLE_CTORS(SpectralBasis, SpectralBasisBase);


#ifdef HAVE_SUNDANCE_STOKHOS
  /** */
  typedef Stokhos::OrthogPolyBasis<int, double> PolyBasis;
  /** */
  typedef Stokhos::OneDOrthogPolyBasis<int, double> PolyBasis1D;

  /** */
  SpectralBasis(const PolyBasis* basis)
    : Handle<SpectralBasisBase>(
      rcp(new StokhosBasisWrapper(rcp(basis)))
      ) 
    {}

  /** */
  SpectralBasis(const PolyBasis1D* basis)
    : Handle<SpectralBasisBase>(rcp(new StokhosBasisWrapper(rcp(basis))))
    {}
#endif

  /** Return the dim of the Spectral Basis */
  int getDim() const {return ptr()->getDim();}

  /** Return the order of the Spectral Basis */
  int getOrder() const {return ptr()->getOrder();}

  /** Return the maximum number of terms */
  int nterms() const {return ptr()->nterms();}
    
  /** Return the basis element stored in the basis array index */
  int getElement(int i) const {return ptr()->getElement(i);}
    
  /** expectation operator */
  double expectation(int i, int j, int k) const 
    {return ptr()->expectation(i,j,k);}

  /** Write to a std::string */
  std::string toString() const {return ptr()->toString();}
};
}

namespace std
{
/** \relates  Sundance::SpectralBasis */
inline ostream& operator<<(std::ostream& os, const Sundance::SpectralBasis& s)
{
  os << s.toString();
  return os;
}
}

#endif
