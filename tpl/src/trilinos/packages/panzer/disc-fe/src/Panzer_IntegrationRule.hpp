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


#ifndef PANZER_INTEGRATION_RULE_HPP
#define PANZER_INTEGRATION_RULE_HPP

#include "Teuchos_ArrayRCP.hpp"
#include "Shards_CellTopology.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_PointRule.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_CubatureControlVolume.hpp"
#include "Intrepid2_CubatureControlVolumeSide.hpp"
#include "Intrepid2_CubatureControlVolumeBoundary.hpp"
#include "Kokkos_DynRankView.hpp"

namespace panzer {

  class CellData;

  /** Derived class for building a point rule based on a finite element integration rule.

      \param[in] cubature_degree Order of the cubature integration.
      \param[in] cell_data Description of the cell.
   */
  class IntegrationRule : public PointRule {
  public:
    
    //! if side = -1 then we use the cell volume integration rule.
    IntegrationRule(int cubature_degree, const panzer::CellData& cell_data);

    IntegrationRule(const panzer::CellData& cell_data, std::string cv_type);

    void setup(int cubature_degree, const panzer::CellData& cell_data);

    void setup_cv(const panzer::CellData& cell_data, std::string cv_type);
  
    //! Returns the order of integration (cubature degree in intrepid lingo)
    int order() const;

    int cubature_degree;

    std::string cv_type;

    //! print information about the integration rule
    virtual void print(std::ostream & os);

    //! Construct an array containing the reference coordinates 
    void referenceCoordinates(Kokkos::DynRankView<double,PHX::Device> & container);
  
  private:

  };

}

#endif
