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


/** \file   Intrepid_OrientationToolsDef.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HDIV_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HDIV_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  namespace Impl {

    template<typename outputViewType,
             typename subcellBasisType,
             typename cellBasisType>
    inline
    void
    OrientationTools::
    getCoeffMatrix_HDIV(outputViewType &output,
                        const subcellBasisType subcellBasis,
                        const cellBasisType cellBasis,
                        const ordinal_type subcellId,
                        const ordinal_type subcellOrt) {
      typedef typename outputViewType::execution_space space_type;
      typedef typename outputViewType::value_type value_type;
      
      // with shards, everything should be computed on host space
      typedef typename
        Kokkos::Impl::is_space<space_type>::host_mirror_space::execution_space host_space_type;

      typedef Kokkos::DynRankView<value_type,host_space_type> DynRankViewHostType;

      ///
      /// Topology
      ///
      
      // populate points on a subcell and map to subcell
      const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      const shards::CellTopology subcellTopo = subcellBasis.getBaseCellTopology();

      const ordinal_type cellDim = cellTopo.getDimension();
      const ordinal_type subcellDim = subcellTopo.getDimension();

      INTREPID2_TEST_FOR_EXCEPTION( subcellDim >= cellDim,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "cellDim must be greater than subcellDim.");

      const auto subcellBaseKey = subcellTopo.getBaseKey();
      //const auto cellBaseKey = cellTopo.getBaseKey();
        
      // if node map has left handed system, orientation should be re-enumerated.
      ordinal_type ort = -1;
      switch (subcellBaseKey) {
      case shards::Line<>::key: {
        if (subcellOrt >= 0 && subcellOrt <  2)
          ort = subcellOrt;
        break;
      }
      case shards::Triangle<>::key: {
        if (subcellOrt >= 0 && subcellOrt <  6) {
          const ordinal_type leftHanded = cellTopo.getNodeMap(2, subcellId, 1) > cellTopo.getNodeMap(2, subcellId, 2);
          const ordinal_type leftOrt[] = { 0, 2, 1, 3, 5, 4 };
          ort = (leftHanded ? leftOrt[subcellOrt] : subcellOrt);
        }
        break;
      }
      case shards::Quadrilateral<>::key: {
        if (subcellOrt >= 0 && subcellOrt <  8) {
          const ordinal_type leftHanded = cellTopo.getNodeMap(2, subcellId, 1) > cellTopo.getNodeMap(2, subcellId, 3);
          const ordinal_type leftOrt[] = { 0, 3, 2, 1, 4, 7, 6, 5 };
          ort = (leftHanded ? leftOrt[subcellOrt] : subcellOrt);
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( subcellBaseKey != shards::Line<>::key ||
                                      subcellBaseKey != shards::Quadrilateral<>::key ||
                                      subcellBaseKey != shards::Triangle<>::key,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                      "subcellBasis must have line, quad, or triangle topology.");
      }
      }
      INTREPID2_TEST_FOR_EXCEPTION( ort == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "Orientation is not properly setup.");
      
      ///
      /// Function space
      ///

      
      {
        const std::string cellBasisName(cellBasis.getName());
        if (cellBasisName.find("HDIV") != std::string::npos) {
          const std::string subcellBasisName(subcellBasis.getName());
          // edge hdiv is hgrad with gauss legendre points
          INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HGRAD") == std::string::npos,
                                        std::logic_error,
                                        ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
                                        "subcellBasis function space is not consistent to cellBasis.");
        }
      }


      ///
      /// Collocation points
      ///
      
      const ordinal_type degree = cellBasis.getDegree();
      
      const ordinal_type numCellBasis = cellBasis.getCardinality();
      const ordinal_type numSubcellBasis = subcellBasis.getCardinality();

      const ordinal_type ordSubcell = cellBasis.getDofOrdinal(subcellDim, subcellId, 0);
      INTREPID2_TEST_FOR_EXCEPTION( ordSubcell == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "Invalid subcellId returns -1 ordSubcell.");

      const ordinal_type ndofSubcell = cellBasis.getDofTag(ordSubcell)(3);

      // reference points on a subcell
      DynRankViewHostType refPtsSubcell;

      switch (subcellBaseKey) {
      case shards::Line<>::key:
      case shards::Triangle<>::key: {
        // hdiv ndof: p-1, the interior points of p+1 with offset 1: p-1.
        const ordinal_type ndofLine = PointTools::getLatticeSize(subcellTopo, degree+1, 1);
        refPtsSubcell = DynRankViewHostType("refPtsSubcell", ndofLine, subcellDim);
        PointTools::getLattice(refPtsSubcell,
                               subcellTopo,
                               degree+1,
                               1, // offset by 1 so the points are located inside
                               POINTTYPE_EQUISPACED);
        INTREPID2_TEST_FOR_EXCEPTION( ndofSubcell != ndofLine,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                      "The number of DOFs in line should be equal to the number of collocation points.");
        break;
      }
      case shards::Quadrilateral<>::key: {
        // hdiv ndof: (p-1)*(p-1), tensor product of interior points of p+1 with offset 1: (p-1)*(p-1).

        // tensor product of lines
        const auto lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
        const ordinal_type ndofLine = PointTools::getLatticeSize(lineTopo, degree+1, 1);
        DynRankViewHostType refPtsLine("refPtsLine", ndofLine, 1);
        PointTools::getLattice(refPtsLine,
                               lineTopo,
                               degree+1,
                               1, // offset by 1 so the points are located inside
                               POINTTYPE_EQUISPACED);

        refPtsSubcell = DynRankViewHostType("refPtsSubcell", ndofLine*ndofLine, subcellDim);
        ordinal_type idx = 0;
        for (ordinal_type j=0;j<ndofLine;++j) { // y
          for (ordinal_type i=0;i<ndofLine;++i,++idx) { // x
            refPtsSubcell(idx, 0) = refPtsLine(i);
            refPtsSubcell(idx, 1) = refPtsLine(j);
          }
        }
        INTREPID2_TEST_FOR_EXCEPTION( idx != ndofSubcell,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                      "counted points on quad is less than ndofSubcell.");
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( subcellBaseKey != shards::Line<>::key ||
                                      subcellBaseKey != shards::Quadrilateral<>::key ||
                                      subcellBaseKey != shards::Triangle<>::key,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                      "subcellBasis must have line, quad, or triangle topology.");
      }
      }
      
      const ordinal_type nptsSubcell = refPtsSubcell.dimension(0);
      
      // modified points with orientation
      DynRankViewHostType ortPtsSubcell("ortPtsSubcell", nptsSubcell, subcellDim);
      Impl::OrientationTools::mapToModifiedReference(ortPtsSubcell,
                                                     refPtsSubcell,
                                                     subcellTopo,
                                                     ort);
      
      // map to reference coordinates
      DynRankViewHostType refPtsCell("refPtsCell", nptsSubcell, cellDim);
      CellTools<host_space_type>::mapToReferenceSubcell(refPtsCell,
                                                        refPtsSubcell,
                                                        subcellDim,
                                                        subcellId,
                                                        cellTopo);

      ///
      /// Basis evaluation on the collocation points
      ///

      // evaluate values on the reference cell
      DynRankViewHostType refValues("refValues", numCellBasis, nptsSubcell, cellDim);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);

      // evaluate values on the modified cell
      DynRankViewHostType outValues("outValues", numSubcellBasis, nptsSubcell);
      subcellBasis.getValues(outValues, ortPtsSubcell);

      auto orient_values = [&](const ordinal_type c) { 
        for (ordinal_type i=0;i<ndofSubcell;++i) {
          const ordinal_type ii = subcellBasis.getDofOrdinal(subcellDim, 0, i);
          for (ordinal_type j=0;j<nptsSubcell;++j) 
            outValues(ii,j) *= c;
        }
      };

      switch (subcellBaseKey) {
      case shards::Line<>::key: {
        // second dimension is dummy
        const ordinal_type c[2] = { 1, -1 };
        orient_values(c[ort]);
        break;
      }
      case shards::Triangle<>::key: {
        const ordinal_type c[6] = {  1,  1,  1,
                                    -1, -1, -1 };
        orient_values(c[ort]);
        break;
      }
      case shards::Quadrilateral<>::key: {
        const ordinal_type c[8] = {  1,  1,  1,  1,
                                    -1, -1, -1, -1 };
        orient_values(c[ort]);
        break;
      }
      }
      
      ///
      /// Restrict vector valued basis functions on the subcell dimensions
      ///
      auto normalize = [&](DynRankViewHostType v) { 
        value_type norm = 0.0;
        const ordinal_type iend = v.dimension(0);
        for (ordinal_type i=0;i<iend;++i)
          norm += v(i)*v(i);
        norm = std::sqrt(norm);
        for (ordinal_type i=0;i<iend;++i)
          v(i) /= norm;
      };
      
      switch (subcellBaseKey) {
      case shards::Line<>::key: {
        DynRankViewHostType edgeNormal("edgeNormal", cellDim);
        CellTools<host_space_type>::getReferenceEdgeTangent(edgeNormal, subcellId, cellTopo);

        normalize(edgeNormal);

        // 90 degree rotation (for 2D)
        if (subcellDim == 2) {
          const auto tan_x = edgeNormal(0);
          const auto tan_y = edgeNormal(1);
          
          edgeNormal(0) = -tan_y;
          edgeNormal(1) =  tan_x;
        }
        // why this does not work ? normal is not defined on 2D side ?
        // DynRankViewHostType edgeNormal("edgeNormal", cellDim);
        // CellTools<host_space_type>::getReferenceSideNormal(edgeNormal, subcellId, cellTopo);

        // normalize(edgeNormal);

        DynRankViewHostType tmpValues("tmpValues", numCellBasis, nptsSubcell);
        for (ordinal_type i=0;i<numCellBasis;++i)
          for (ordinal_type j=0;j<nptsSubcell;++j)
            for (ordinal_type k=0;k<cellDim;++k)
              tmpValues(i,j) += refValues(i,j,k)*edgeNormal(k);
        refValues = tmpValues;
        break;
      }
      case shards::Quadrilateral<>::key: {
        DynRankViewHostType faceNormal("faceNormal", cellDim);
        CellTools<host_space_type>::getReferenceFaceNormal(faceNormal, subcellId, cellTopo);

        normalize(faceNormal);

        DynRankViewHostType tmpValues("tmpValues", numCellBasis, nptsSubcell);
        for (ordinal_type i=0;i<numCellBasis;++i)
          for (ordinal_type j=0;j<nptsSubcell;++j)
            for (ordinal_type k=0;k<cellDim;++k) 
              tmpValues(i,j) += refValues(i,j,k)*faceNormal(k);
        refValues = tmpValues;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, "Should not come here" );        
      }
      }
      
      ///
      /// Construct collocation matrix and solve problems
      ///

      // construct collocation matrix; using lapack, it should be left layout
      Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space_type>
        refMat("refMat", nptsSubcell, ndofSubcell),
        ortMat("ortMat", nptsSubcell, ndofSubcell),
        pivVec("pivVec", nptsSubcell, 1);

      for (ordinal_type i=0;i<ndofSubcell;++i) {
        const ordinal_type iref = cellBasis.getDofOrdinal(subcellDim, subcellId, i);
        const ordinal_type iout = subcellBasis.getDofOrdinal(subcellDim, 0, i);

        for (ordinal_type j=0;j<nptsSubcell;++j) {
          refMat(j,i) = refValues(iref,j);
          ortMat(j,i) = outValues(iout,j);
        }
      }
      
      // solve the system
      {
        Teuchos::LAPACK<ordinal_type,value_type> lapack;
        ordinal_type info = 0;

        lapack.GESV(nptsSubcell, ndofSubcell,
                    refMat.data(),
                    refMat.stride_1(),
                    (ordinal_type*)pivVec.data(),
                    ortMat.data(),
                    ortMat.stride_1(),
                    &info);

        if (info) {
          std::stringstream ss;
          ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
             << "LAPACK return with error code: "
             << info;
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
        }
        
        // clean up numerical noise
        {
          const double eps = threshold();
          const ordinal_type
            iend = ortMat.dimension(0),
            jend = ortMat.dimension(1);
          for (ordinal_type i=0;i<iend;++i)
            for (ordinal_type j=0;j<jend;++j)
              if (std::abs(ortMat(i,j)) < eps) ortMat(i,j) = 0;
        }
      }

      {
        // move the data to original device memory
        const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofSubcell);
        Kokkos::deep_copy(Kokkos::subview(output, range, range),
                          Kokkos::subview(ortMat, range, range));
      }
    }
  }

}
#endif
