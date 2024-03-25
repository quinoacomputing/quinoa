// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Extract field output for inciter
  \details   Extract field output for inciter.
*/
// *****************************************************************************

#include "FieldOutput.hpp"
#include "ContainerUtil.hpp"
#include "Vector.hpp"
#include "Integrate/Basis.hpp"

namespace inciter {

std::vector< std::string >
numericFieldNames( tk::Centering c, char /*depvar*/ )
// *****************************************************************************
// Collect field output names from numerical solution based on user input
//! \param[in] c Extract variable names only with this centering
// //! \param[in] depvar Consider this depvar (mesh) only, ignore if 0
//! \return Output field names requested by user
// *****************************************************************************
{
  std::vector< std::string > f;
  for (const auto& v : g_inputdeck.get< newtag::field_output, newtag::outvar >()) {
    if (v.centering == c && !v.analytic()) {
      std::stringstream s;
      s << v.name;
      f.push_back( s.str() );
    }
  }

  return f;
}

std::vector< std::vector< tk::real > >
numericFieldOutput( const tk::Fields& U,
                    tk::Centering c,
                    const tk::Fields& P,
                    char /*depvar*/ )
// *****************************************************************************
// Collect field output from numerical solution based on user input
//! \param[in] U Solution data to extract from
//! \param[in] c Extract variables only with this centering
//! \param[in] P Optional primitive variable solution data to extract from
// //! \param[in] depvar Consider this depvar (mesh) only, ignore if 0
//! \return Output fields requested by user
// *****************************************************************************
{
  // will not use P if empty
  const auto& p = P.empty() ? U : P;

  //auto rdof =
  //  c == tk::Centering::NODE ? 1 : g_inputdeck.get< newtag::rdof >();
  std::size_t rdof = 1;

  std::vector< std::vector< tk::real > > f;
  for (const auto& v : g_inputdeck.get< newtag::field_output, newtag::outvar >()) {
    if (v.centering == c) {
      const auto& F = v.primitive() ? p : U;
      if (v.name.empty()) {        // depvar-based direct access
        f.push_back( F.extract_comp( v.field*rdof ) );
      } else if (!v.analytic()) {  // human-readable non-analytic via custom fn
        Assert( v.getvar, "getvar() not configured for " + v.name );
        f.push_back( v.getvar( F, rdof ) );
      }
    }
  }

  return f;
}

void
evalSolution(
  const Discretization& D,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, std::size_t >& addedTets,
  const std::vector< std::size_t >& ndofel,
  const tk::Fields& U,
  const tk::Fields& P,
  tk::Fields& uElemfields,
  tk::Fields& pElemfields,
  tk::Fields& uNodefields,
  tk::Fields& pNodefields )
// *****************************************************************************
// Evaluate solution on incoming (a potentially refined) mesh
//! \param[in] D Discretization base class to read from
//! \param[in] inpoel Incoming (potentially refined field-output) mesh
//!   connectivity
//! \param[in] coord Incoming (potentially refined Field-output) mesh node
//!   coordinates
//! \param[in] addedTets Field-output mesh cells and their parents (local ids)
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] U Solution vector
//! \param[in] P Vector of primitives
//! \param[in,out] uElemfields Solution elem output fields
//! \param[in,out] pElemfields Primitive elem output fields
//! \param[in,out] uNodefields Solution nodal output fields
//! \param[in,out] pNodefields Primitive nodal output fields
//! \details This function evaluates the solution on the incoming mesh, and
//!   stores it in uElemfields, pElemfields, uNodefields, and pNodefields
//!   appropriately. The incoming mesh can be refined but can also be just the
//!   mesh the numerical solution is computed on.
//! \note If the incoming mesh is refined (for field output) compared to the
//!   mesh the numerical solution is computed on, the solution is evaluated in
//!   cells as wells as in nodes. If the solution is not refined, the solution
//!   is evaluated in nodes.
// *****************************************************************************
{
  using tk::dot;
  using tk::real;

  const auto nelem = inpoel.size()/4;
  const auto rdof = g_inputdeck.get< newtag::rdof >();
  const auto uncomp = U.nprop() / rdof;
  const auto pncomp = P.nprop() / rdof;

  // If mesh is not refined for field output, cut off ghosts from element
  // solution. (No need to output ghosts and writer would error.) If mesh is
  // refined for field output, resize element solution fields to refined mesh.
  uElemfields.resize( nelem );
  pElemfields.resize( nelem );

  auto npoin = coord[0].size();
  uNodefields.resize( npoin );
  pNodefields.resize( npoin );
  uNodefields.fill(0.0);
  pNodefields.fill(0.0);

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Assign values to element-fields
  for (std::size_t e=0; e<U.nunk(); ++e) {
    if (e < nelem) {
      for (std::size_t i=0; i<uncomp; ++i) {
        uElemfields(e,i) = U(e,rdof*i);
      }
      for (std::size_t i=0; i<pncomp; ++i) {
        pElemfields(e,i) = P(e,rdof*i);
      }
    }
  }

  // If mesh is not refined for output, evaluate solution in nodes
  if (addedTets.empty()) {

    for (std::size_t e=0; e<nelem; ++e) {
      std::size_t dofe(1);
      if (!ndofel.empty()) {
        dofe = ndofel[e];
      }
      auto e4 = e*4;
      // Extract element node coordinates
      std::array< std::array< real, 3>, 4 > ce{{
        {{ x[inpoel[e4  ]], y[inpoel[e4  ]], z[inpoel[e4  ]] }},
        {{ x[inpoel[e4+1]], y[inpoel[e4+1]], z[inpoel[e4+1]] }},
        {{ x[inpoel[e4+2]], y[inpoel[e4+2]], z[inpoel[e4+2]] }},
        {{ x[inpoel[e4+3]], y[inpoel[e4+3]], z[inpoel[e4+3]] }} }};
      // Compute inverse Jacobian
      auto J = tk::inverseJacobian( ce[0], ce[1], ce[2], ce[3] );
      // Evaluate solution in child nodes
      for (std::size_t j=0; j<4; ++j) {
        std::array< real, 3 >
           h{{ce[j][0]-ce[0][0], ce[j][1]-ce[0][1], ce[j][2]-ce[0][2] }};
        auto Bn = tk::eval_basis( dofe,
                                  dot(J[0],h), dot(J[1],h), dot(J[2],h) );
        auto u = eval_state( uncomp, rdof, dofe, e, U, Bn );
        auto p = eval_state( pncomp, rdof, dofe, e, P, Bn );
        // Assign child node solution
        for (std::size_t i=0; i<uncomp; ++i) uNodefields(inpoel[e4+j],i) += u[i];
        for (std::size_t i=0; i<pncomp; ++i) pNodefields(inpoel[e4+j],i) += p[i];
      }
    }

  // If mesh is refed for output, evaluate solution in elements and nodes of
  // refined mesh
  } else {

    const auto& pinpoel = D.Inpoel();  // unrefined (parent) mesh

    for ([[maybe_unused]] const auto& [child,parent] : addedTets) {
      Assert( child < nelem, "Indexing out of new solution vector" );
      Assert( parent < pinpoel.size()/4,
              "Indexing out of old solution vector" );
    }

    for (const auto& [child,parent] : addedTets) {
      std::size_t dofe(1);
      if (!ndofel.empty()) {
        dofe = ndofel[parent];
      }
      // Extract parent element's node coordinates
      auto p4 = 4*parent;
      std::array< std::array< real, 3>, 4 > cp{{
        {{ x[pinpoel[p4  ]], y[pinpoel[p4  ]], z[pinpoel[p4  ]] }},
        {{ x[pinpoel[p4+1]], y[pinpoel[p4+1]], z[pinpoel[p4+1]] }},
        {{ x[pinpoel[p4+2]], y[pinpoel[p4+2]], z[pinpoel[p4+2]] }},
        {{ x[pinpoel[p4+3]], y[pinpoel[p4+3]], z[pinpoel[p4+3]] }} }};
      // Evaluate inverse Jacobian of the parent
      auto Jp = tk::inverseJacobian( cp[0], cp[1], cp[2], cp[3] );
      // Compute child cell centroid
      auto c4 = 4*child;
      auto cx = (x[inpoel[c4  ]] + x[inpoel[c4+1]] +
                 x[inpoel[c4+2]] + x[inpoel[c4+3]]) / 4.0;
      auto cy = (y[inpoel[c4  ]] + y[inpoel[c4+1]] +
                 y[inpoel[c4+2]] + y[inpoel[c4+3]]) / 4.0;
      auto cz = (z[inpoel[c4  ]] + z[inpoel[c4+1]] +
                 z[inpoel[c4+2]] + z[inpoel[c4+3]]) / 4.0;
      // Compute solution in child centroid
      std::array< real, 3 > h{{cx-cp[0][0], cy-cp[0][1], cz-cp[0][2] }};
      auto B = tk::eval_basis( dofe, dot(Jp[0],h), dot(Jp[1],h), dot(Jp[2],h) );
      auto u = eval_state( uncomp, rdof, dofe, parent, U, B );
      auto p = eval_state( pncomp, rdof, dofe, parent, P, B );
      // Assign cell center solution from parent to child
      for (std::size_t i=0; i<uncomp; ++i) uElemfields(child,i) = u[i];
      for (std::size_t i=0; i<pncomp; ++i) pElemfields(child,i) = p[i];
      // Extract child element's node coordinates
      std::array< std::array< real, 3>, 4 > cc{{
        {{ x[inpoel[c4  ]], y[inpoel[c4  ]], z[inpoel[c4  ]] }},
        {{ x[inpoel[c4+1]], y[inpoel[c4+1]], z[inpoel[c4+1]] }},
        {{ x[inpoel[c4+2]], y[inpoel[c4+2]], z[inpoel[c4+2]] }},
        {{ x[inpoel[c4+3]], y[inpoel[c4+3]], z[inpoel[c4+3]] }} }};
      // Evaluate solution in child nodes
      for (std::size_t j=0; j<4; ++j) {
        std::array< real, 3 >
           hn{{cc[j][0]-cp[0][0], cc[j][1]-cp[0][1], cc[j][2]-cp[0][2] }};
        auto Bn = tk::eval_basis( dofe,
                                  dot(Jp[0],hn), dot(Jp[1],hn), dot(Jp[2],hn) );
        auto cnu = eval_state(uncomp, rdof, dofe, parent, U, Bn);
        auto cnp = eval_state(pncomp, rdof, dofe, parent, P, Bn);
        // Assign child node solution
        for (std::size_t i=0; i<uncomp; ++i) uNodefields(inpoel[c4+j],i) += cnu[i];
        for (std::size_t i=0; i<pncomp; ++i) pNodefields(inpoel[c4+j],i) += cnp[i];
      }
    }
  }
}

} // inciter::
