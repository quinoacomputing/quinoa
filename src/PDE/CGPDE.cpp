// *****************************************************************************
/*!
  \file      src/PDE/CGPDE.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions common to ALECG
  \details   Functions common to ALECG.
*/
// *****************************************************************************

#include <array>
#include <vector>
#include <unordered_map>

#include "Vector.hpp"
#include "DerivedData.hpp"
#include "Exception.hpp"
#include "Around.hpp"
#include "Fields.hpp"
#include "CGPDE.hpp"
#include "FunctionPrototypes.hpp"

namespace inciter {
namespace cg {

std::vector< tk::real >
solinc( tk::ncomp_t ncomp,
        const std::vector< EOS >& mat_blk, tk::real x, tk::real y,
        tk::real z, tk::real t, tk::real dt, tk::InitializeFn solution )
// *****************************************************************************
// Evaluate the increment from t to t+dt of the analytical solution at (x,y,z)
// for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution increment starting from
//! \param[in] dt Time increment at which evaluate the solution increment to
//! \param[in] solution Function used to evaluate the solution
//! \return Increment in values of all components evaluated at (x,y,z,t+dt)
// *****************************************************************************
{
  auto st1 = solution( ncomp, mat_blk, x, y, z, t );
  auto st2 = solution( ncomp, mat_blk, x, y, z, t+dt );

  std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                  []( tk::real s, tk::real& d ){ return d -= s; } );

  return st2;
}

std::unordered_map< int,
  std::unordered_map< std::size_t, std::array< tk::real, 4 > > >
bnorm(
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< int, std::unordered_set< std::size_t > >& bcnodes )
// *****************************************************************************
//! Compute boundary point normals
//! \param[in] bface Boundary-faces mapped to side sets used in the input file
//! \param[in] triinpoel Boundary-face connectivity where BCs set (local ids)
//! \param[in] coord Mesh node coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bcnodes Local node ids associated to side set ids at which BCs
//!   are set that require normals
//! \return Face normals in boundary points, Inner key: local node id, value:
//!   unit normal and inverse distance square between face centroids and points,
//!   outer key: side set id
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Lambda to compute the inverse distance squared between boundary face
  // centroid and boundary point. Here p is the global node id and g is the
  // geometry of the boundary face, see tk::geoFaceTri().
  auto invdistsq = [&]( const tk::Fields& g, std::size_t p ){
    return 1.0 / ( (g(0,4) - x[p])*(g(0,4) - x[p]) +
                   (g(0,5) - y[p])*(g(0,5) - y[p]) +
                   (g(0,6) - z[p])*(g(0,6) - z[p]) );
  };

  // Compute boundary point normals on all side sets summing inverse distance
  // weighted face normals to points. This is only a partial sum at shared
  // boundary points in parallel. Inner key: global node id, value: normals and
  // inverse distance square, outer key, side set id.
  std::unordered_map< int,
    std::unordered_map< std::size_t, std::array< tk::real, 4 > > > norm;
  for (const auto& [ setid, faceids ] : bface) { // for all side sets
    for (auto f : faceids) { // for all side set triangles
      tk::UnsMesh::Face
        face{ triinpoel[f*3+0], triinpoel[f*3+1], triinpoel[f*3+2] };
      std::array< tk::real, 3 > fx{ x[face[0]], x[face[1]], x[face[2]] };
      std::array< tk::real, 3 > fy{ y[face[0]], y[face[1]], y[face[2]] };
      std::array< tk::real, 3 > fz{ z[face[0]], z[face[1]], z[face[2]] };
      auto g = tk::geoFaceTri( fx, fy, fz );
      for (auto p : face) {  // for all 3 nodes of a boundary triangle face
        for (const auto& [s,nodes] : bcnodes) {  // for all bnd nodes w/ normals
          if (setid == s) {  // only contribute to side set we operate on
            auto i = nodes.find(p);
            if (i != end(nodes)) {        // only if user set bc on node
              tk::real r = invdistsq(g,p);
              auto& n = norm[s][gid[p]];  // associate set id and global node id
              n[0] += r*g(0,1);
              n[1] += r*g(0,2);
              n[2] += r*g(0,3);
              n[3] += r;
            }
          }
        }
      }
    }
  }

  return norm;
}

} // cg::
} // inciter::
