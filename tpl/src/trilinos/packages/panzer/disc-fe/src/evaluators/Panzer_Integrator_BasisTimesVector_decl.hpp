// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_Integerator_BasisTimesVector_decl_hpp__
#define __Panzer_Integerator_BasisTimesVector_decl_hpp__

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
PANZER_EVALUATOR_CLASS(Integrator_BasisTimesVector)
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
    
  PHX::MDField<ScalarT,Cell,IP,Dim> vectorField;

  PHX::MDField<ScalarT,Cell,BASIS> dof_orientation;

  std::vector<PHX::MDField<ScalarT,Cell,IP> > field_multipliers;
  Kokkos::View<Kokkos::View<ScalarT** >* > kokkos_field_multipliers;

  std::size_t basis_card;

  std::size_t num_qp;

  std::size_t num_dim;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

public:

  typedef PHX::Device::scratch_memory_space shmem_space ;

  template<int NUM_FIELD_MULT>
  struct FieldMultTag{};

  template<int NUM_FIELD_MULT>
  void operator()(const FieldMultTag<NUM_FIELD_MULT> &, const std::size_t &cell) const;
private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  PHX::MDField<double,Cell,BASIS,IP,Dim> weighted_basis_vector;

PANZER_EVALUATOR_CLASS_END

}

#endif
