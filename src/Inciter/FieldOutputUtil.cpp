// *****************************************************************************
/*!
  \file      src/Inciter/FieldOutputUtil.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field output utility functionality
  \details   Field output utility functions.
*/
// *****************************************************************************

#include "FieldOutputUtil.hpp"
#include "Vector.hpp"
#include "Around.hpp"
#include "Integrate/Basis.hpp"

std::tuple< tk::Fields, tk::Fields >
tk::nodeEval( std::size_t offset,
              std::size_t ndof,
              std::size_t rdof,
              const tk::UnsMesh::Coords& coord,
              const std::vector< std::size_t >& inpoel,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& esup,
              const Fields& U,
              const Fields& P )
// *****************************************************************************
//  Evaluate solution in nodes
//! \param[in] offset Index for equation systems
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] coord Node coordinates
//! \param[in] inpoel Mesh connectivity
//! \param[in] esup Elements surrounding points
//! \param[in] U Vector of cell-averaged unknowns
//! \param[in] P Vector of cell-averaged primitive quantities
//! \return Vector of unknowns at nodes, vector of primitive quantities at nodes
// *****************************************************************************
{
  using tk::dot;

  auto npoin = coord[0].size();

  tk::Fields Un( npoin, U.nprop()/rdof );
  tk::Fields Pn( npoin, P.nprop()/rdof );

  Un.fill(0.0);
  Pn.fill(0.0);

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t p=0; p<npoin; ++p) {
    for (auto e : tk::Around(esup,p)) {
      // Extract element coordinates
      auto e4 = 4*e;
      std::array< std::array< tk::real, 3>, 4 > ce{{
        {{ x[inpoel[e4  ]], y[inpoel[e4  ]], z[inpoel[e4  ]] }},
        {{ x[inpoel[e4+1]], y[inpoel[e4+1]], z[inpoel[e4+1]] }},
        {{ x[inpoel[e4+2]], y[inpoel[e4+2]], z[inpoel[e4+2]] }},
        {{ x[inpoel[e4+3]], y[inpoel[e4+3]], z[inpoel[e4+3]] }} }};
      // Evaluate inverse Jacobian of the element
      auto J = tk::inverseJacobian( ce[0], ce[1], ce[2], ce[3] );
      // Evaluate in and sum solution to p (will average later)
      std::array<tk::real,3> h{ x[p]-ce[0][0], y[p]-ce[0][1], z[p]-ce[0][2] };
      auto B = tk::eval_basis( ndof, dot(J[0],h), dot(J[1],h), dot(J[2],h) );
      auto chu = eval_state( Un.nprop(), offset, rdof, ndof, e, U, B );
      for (std::size_t i=0; i<Un.nprop(); ++i) Un(p,i,offset) += chu[i];
      auto chp = eval_state( Pn.nprop(), offset, rdof, ndof, e, P, B );
      for (std::size_t i=0; i<Pn.nprop(); ++i) Pn(p,i,offset) += chp[i];
    }
  }

  return { Un, Pn };
}
