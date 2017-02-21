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

#ifndef PANZER_POINT_VALUES2_HPP
#define PANZER_POINT_VALUES2_HPP

#include "PanzerDiscFE_config.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_Dimension.hpp"

#include "Teuchos_RCP.hpp"

#include "Intrepid2_CellTools.hpp"

namespace panzer {

  template <typename Scalar,
            template <typename DataT,
               typename Tag0, typename Tag1, typename Tag2,
               typename Tag3, typename Tag4, typename Tag5,
               typename Tag6, typename Tag7> class Array >
  struct PointValues2 {
    typedef typename ArrayTraits<Scalar, Array<Scalar,void,void,void,void,void,void,void,void> >::size_type size_type;
    
    //! Sizes/allocates memory for arrays
    template <typename ArrayFactory>
    void setupArrays(const Teuchos::RCP<const panzer::PointRule>& pr,const ArrayFactory & af);

    template <typename NodeCoordinateArray,typename PointCoordinateArray>
    inline void evaluateValues(const NodeCoordinateArray & node_coordinates,const PointCoordinateArray & point_coordinates);

    template <typename CoordinateArray>
    void copyNodeCoords(const CoordinateArray& in_node_coords);

    template <typename CoordinateArray>
    void copyPointCoords(const CoordinateArray& in_point_coords);

    Array<Scalar,IP,Dim,void,void,void,void,void,void> coords_ref;      // <IP,Dim>
    Array<Scalar,Cell,NODE,Dim,void,void,void,void,void> node_coordinates; // <Cell,NODE,Dim>
    Array<Scalar,Cell,IP,Dim,Dim,void,void,void,void> jac;              // <Cell,IP,Dim,Dim>
    Array<Scalar,Cell,IP,Dim,Dim,void,void,void,void> jac_inv;          // <Cell,IP,Dim,Dim>
    Array<Scalar,Cell,IP,void,void,void,void,void,void> jac_det;        // <Cell,IP>

    // cell points
    Array<Scalar,Cell,IP,Dim,void,void,void,void,void> point_coords;    // <Cell,IP,Dim>

    Teuchos::RCP<const panzer::PointRule> point_rule;
  };

  template <typename Scalar,
            template <typename DataT,
               typename Tag0, typename Tag1, typename Tag2,
               typename Tag3, typename Tag4, typename Tag5,
               typename Tag6, typename Tag7> class Array >
  template <typename NodeCoordinateArray,typename PointCoordinateArray>
  void PointValues2<Scalar,Array>::
  evaluateValues(const NodeCoordinateArray& in_node_coords,
                 const PointCoordinateArray & in_point_coords)
  {
    if (point_rule->isSide()) {
       TEUCHOS_ASSERT(false); // not implemented!!!!
    }
    
    copyPointCoords(in_point_coords);
    copyNodeCoords(in_node_coords);
    
    Intrepid2::CellTools<Scalar> cell_tools;
    
    cell_tools.setJacobian(jac, coords_ref, node_coordinates,*(point_rule->topology));
    cell_tools.setJacobianInv(jac_inv, jac);
    cell_tools.setJacobianDet(jac_det, jac);
    
    // IP coordinates
    cell_tools.mapToPhysicalFrame(point_coords, coords_ref, node_coordinates, *(point_rule->topology));
  }

} // namespace panzer

#include "Panzer_PointValues2_impl.hpp"

#endif
