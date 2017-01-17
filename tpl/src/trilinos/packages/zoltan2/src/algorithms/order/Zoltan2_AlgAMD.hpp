// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_AlgAMD.hpp
    \brief The AMD ordering algorithm uses SuiteSparse.
*/

#ifndef _ZOLTAN2_ALGAMD_HPP_
#define _ZOLTAN2_ALGAMD_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


#ifdef HAVE_ZOLTAN2_AMD
#include "amd.h"
#endif

namespace Zoltan2{

#ifdef HAVE_ZOLTAN2_AMD
template <typename Ordinal>
class AMDTraits
{
    public:
    Ordinal order(Ordinal n, const Ordinal *Ap, const Ordinal *Ai,
                Ordinal *perm, double *control, double *info);
};

template <>
class AMDTraits<int>
{
    public:
    int order(int n, const int *Ap, const int *Ai, int *perm,
                double *control, double *info)
    {
        return (amd_order(n, Ap, Ai, perm, control, info));
    }
};

template <>
class AMDTraits<long>
{
    public:
    long order(long n, const long *Ap, const long *Ai, long *perm,
                double *control, double *info)
    {
        return (amd_l_order(n, Ap, Ai, perm, control, info));
    }
};
#endif

}


namespace Zoltan2{

template <typename Adapter>
class AlgAMD : public Algorithm<Adapter>
{
    private:

    const RCP<GraphModel<Adapter> > model;
    const RCP<Teuchos::ParameterList> pl;
    const RCP<const Teuchos::Comm<int> > comm;
      
    public:

    AlgAMD(
      const RCP<GraphModel<Adapter> > &model__,
      const RCP<Teuchos::ParameterList> &pl__,
      const RCP<const Teuchos::Comm<int> > &comm__
    ) : model(model__), pl(pl__), comm(comm__)
    { }

    int order(const RCP<OrderingSolution<typename Adapter::lno_t,
                                         typename Adapter::gno_t> > &solution)
    {
#ifndef HAVE_ZOLTAN2_AMD
  throw std::runtime_error(
        "BUILD ERROR: AMD requested but not compiled into Zoltan2.\n"
        "Please set CMake flag Zoltan2_ENABLE_AMD:BOOL=ON.");
#else
      typedef typename Adapter::gno_t gno_t;
      typedef typename Adapter::lno_t lno_t;
      typedef typename Adapter::scalar_t scalar_t;

      int ierr= 0;

      const size_t nVtx = model->getLocalNumVertices();

      //cout << "Local num vertices" << nVtx << endl;
      ArrayView<const gno_t> edgeIds;
      ArrayView<const lno_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > wgts;

      // wgts are ignored in AMD
      model->getEdgeList(edgeIds, offsets, wgts);

      AMDTraits<lno_t> AMDobj;
      double Control[AMD_CONTROL];
      double Info[AMD_INFO];

      amd_defaults(Control);
      amd_control(Control);

      lno_t *perm;
      perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());

      lno_t result = AMDobj.order(nVtx, offsets.getRawPtr(),
                             edgeIds.getRawPtr(), perm, Control, Info);

      if (result != AMD_OK && result != AMD_OK_BUT_JUMBLED)
          ierr = -1;

      solution->setHavePerm(true);
      return ierr;
#endif
    }
};

}



#endif
