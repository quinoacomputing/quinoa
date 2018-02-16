// *****************************************************************************
/*!
  \file      src/Mesh/Gradients.C
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Functions computing gradients on unstructured meshes for tetrahedra
  \details   Functions computing gradients using linear finite element shape
             functions on unstructured meshes for tetrahedra.
*/
// *****************************************************************************

#include "Exception.h"
#include "Gradients.h"
#include "Vector.h"

namespace tk {

std::array< tk::real, 3 >
nodegrad( std::size_t node,
          const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup,
          const tk::Fields& U,
          ncomp_t c )
// *****************************************************************************
//  Compute gradient at a mesh node
//! \param[in] node Node id at which to compute gradient
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esup Linked lists storing elements surrounding points, see
//!    tk::genEsup()
//! \param[in] U Field vector whose component gradient to compute
//! \param[in] c Scalar component to compute gradient of
//! \return Gradient of U(c) at mesh node
// *****************************************************************************
{
  Assert( c < U.nprop(), "Indexing out of field data" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // storage for gradient and volume at the mesh node
  std::array< tk::real, 3 > g{{ 0.0, 0.0, 0.0 }};
  tk::real vol = 0.0;

  // loop over cells surrounding mesh node
  for (auto k=esup.second[node]+1; k<=esup.second[node+1]; ++k) {

     // access element id
     auto e = esup.first[k];

     // access node IDs
     const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                            inpoel[e*4+2], inpoel[e*4+3] }};

     // compute element Jacobi determinant
     const std::array< tk::real, 3 >
       ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
       ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
       da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
     const auto J = tk::triple( ba, ca, da );        // J = 6V
     Assert( J > 0, "Element Jacobian non-positive" );

     // shape function derivatives, nnode*ndim [4][3]
     std::array< std::array< tk::real, 3 >, 4 > grad;
     grad[1] = tk::crossdiv( ca, da, J );
     grad[2] = tk::crossdiv( da, ba, J );
     grad[3] = tk::crossdiv( ba, ca, J );
     for (std::size_t i=0; i<3; ++i)
       grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

     // access field data for scalar component c at nodes of element
     auto u = U.extract( c, 0, N );

     // compute nodal volume: every element contributes their volume / 4
     vol += 5.0*J/120.0;

     // compute gradient over element weighed by cell volume / 4 and sum to node
     for (std::size_t j=0; j<3; ++j) {
       for (std::size_t i=0; i<4; ++i) {
         g[j] += grad[i][j] * u[i] * 5.0*J/120.0;
       }
     }
   }

   // divide components of nodal gradient by nodal volume
   for (std::size_t j=0; j<3; ++j) g[j] /= vol;

   return g;
}

std::array< tk::real, 3 >
edgegrad( std::size_t edge,
          const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esued,
          const tk::Fields& U,
          ncomp_t c )
// *****************************************************************************
//  Compute gradient at a mesh edge
//! \param[in] edgeEdge id at which to compute gradient
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esued Linked lists storing elements surrounding edges, see
//!    tk::genEsued()
//! \param[in] U Field vector whose component gradient to compute
//! \param[in] c Scalar component to compute gradient of
//! \return Gradient of U(c) at mesh edge
// *****************************************************************************
{
  Assert( c < U.nprop(), "Indexing out of field data" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // storage for gradient and volume at the mesh edge
  std::array< tk::real, 3 > g{{ 0.0, 0.0, 0.0 }};
  tk::real vol = 0.0;

  // loop over cells surrounding mesh edge
  for (auto k=esued.second[edge]+1; k<=esued.second[edge+1]; ++k) {

     // access element id
     auto e = esued.first[k];

     // access node IDs
     const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                            inpoel[e*4+2], inpoel[e*4+3] }};

     // compute element Jacobi determinant
     const std::array< tk::real, 3 >
       ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
       ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
       da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
     const auto J = tk::triple( ba, ca, da );        // J = 6V
     Assert( J > 0, "Element Jacobian non-positive" );

     // shape function derivatives, nnode*ndim [4][3]
     std::array< std::array< tk::real, 3 >, 4 > grad;
     grad[1] = tk::crossdiv( ca, da, J );
     grad[2] = tk::crossdiv( da, ba, J );
     grad[3] = tk::crossdiv( ba, ca, J );
     for (std::size_t i=0; i<3; ++i)
       grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

     // access field data for scalar component c at nodes of element
     auto u = U.extract( c, 0, N );

     // compute edge volume: every element contributes their volume / 6
     vol += J/36.0;

     // compute gradient over element weighed by cell volume / 6 and sum to edge
     for (std::size_t j=0; j<3; ++j) {
       for (std::size_t i=0; i<4; ++i) {
         g[j] += grad[i][j] * u[i] * J/36.0;
       }
     }
   }

   // divide components of gradient by edge volume
   for (std::size_t j=0; j<3; ++j) g[j] /= vol;

   return g;
}

} // tk::
