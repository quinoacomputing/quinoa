// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_TRI_Cn_FEM.hpp
    \brief  Header file for the Intrepid2::HGRAD_TRI_Cn_FEM_ORTH class.
    \author Created by Robert Kirby
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_HPP__
#define __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_HPP__

#include "Intrepid2_Basis.hpp"

namespace Intrepid2 {
  
  /** \class  Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH
      \brief  Implementation of the default H(grad)-compatible orthogonal basis (Dubiner) of
      arbitrary degree on triangle.
      
      \remarks
      
      \li   All degrees of freedom are considered to be internal (ie not assembled)
  */
  
  namespace Impl {

    template<ordinal_type n>
    struct OrthPolynomial {
      template<ordinal_type maxOrder,
               ordinal_type maxNumPts,
               typename outputViewType,
               typename inputViewType> 
      KOKKOS_INLINE_FUNCTION
      static void 
      generate( /**/  outputViewType output,
                const inputViewType input,
                const ordinal_type p );
    };

    class Basis_HGRAD_TRI_Cn_FEM_ORTH {
    public:

      template<EOperator opType>
      struct Serial {
        template<typename outputViewType,
                 typename inputViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues( /**/  outputViewType output,
                   const inputViewType  input,
                   const ordinal_type   order,
                   const ordinal_type   opDn = 0 );
      };
      
      template<typename ExecSpaceType, ordinal_type numPtsPerEval,
               typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      static void
      getValues(  /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                  const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                  const ordinal_type order,
                  const EOperator operatorType );
      
    };

  }

  template<typename ExecSpaceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_TRI_Cn_FEM_ORTH
    : public Basis<ExecSpaceType,outputValueType,pointValueType> {
  public:
    typedef double value_type;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_1d_host ordinal_type_array_1d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_2d_host ordinal_type_array_2d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_3d_host ordinal_type_array_3d_host;
    
    /** \brief  Constructor.
     */
    Basis_HGRAD_TRI_Cn_FEM_ORTH( const ordinal_type order );

    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::outputViewType outputViewType;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::pointViewType  pointViewType;
    
    virtual
    void
    getValues( /**/  outputViewType outputValues,
               const pointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const {
#ifdef HAVE_INTREPID2_DEBUG
      Intrepid2::getValues_HGRAD_Args(outputValues,
                                      inputPoints,
                                      operatorType,
                                      this->getBaseCellTopology(),
                                      this->getCardinality() );
#endif
      constexpr ordinal_type numPtsPerEval = 1;
      Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
        getValues<ExecSpaceType,numPtsPerEval>( outputValues, 
                                                inputPoints, 
                                                this->getDegree(),
                                                operatorType );
    }
  };



}// namespace Intrepid2

#include "Intrepid2_HGRAD_TRI_Cn_FEM_ORTHDef.hpp"

#endif

