// $Id$
// $Source$
// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Christopher W. Miller (cmiller@math.umd.edu)
//
// ***********************************************************************
// @HEADER

#include "sandia_sgmga.hpp"
#include "sandia_rules.hpp"
#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::AnisoSparseGridQuadrature<ordinal_type, value_type>::
AnisoSparseGridQuadrature(
  const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis,
  ordinal_type sparse_grid_level, value_type dim_weights[],
  value_type duplicate_tol,
  ordinal_type growth_rate) :
  coordinate_bases(product_basis->getCoordinateBases())
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::AnisoSparseGridQuadrature -- Quad Grid Generation");
#endif
  ordinal_type d = product_basis->dimension();
  ordinal_type p = product_basis->order();
  ordinal_type sz = product_basis->size();
  ordinal_type level = sparse_grid_level;

  // Mike's heuristic formula for computing the level
  if (level == 0) {
    level = static_cast<ordinal_type>(std::ceil(0.5*(p+d-1)));
    if (level < d)
      level = p;
  }

  std::cout << "Sparse grid level = " << level << std::endl;

  // Compute quad points, weights, values
   Teuchos::Array<typename OneDOrthogPolyBasis<ordinal_type,value_type>::LevelToOrderFnPtr> growth_rules(d);

  Teuchos::Array< void (*) ( int order, int dim, double x[] ) > compute1DPoints(d);
  Teuchos::Array< void (*) ( int order, int dim, double w[] ) > compute1DWeights(d);
  for (ordinal_type i=0; i<d; i++) {
    compute1DPoints[i] = &(getMyPoints);
    compute1DWeights[i] = &(getMyWeights);
    growth_rules[i] = coordinate_bases[i]->getSparseGridGrowthRule();
  }

  // Set the static sparse grid quadrature pointer to this
  // (this will cause a conflict if another sparse grid quadrature object
  // is trying to access the VPISparseGrid library, but that's all we can
  // do at this point).
  sgq = this;

  int num_total_pts =
    webbur::sandia_sgmga_size_total(d,&dim_weights[0], level, growth_rate,
				    &growth_rules[0]);

  ordinal_type ntot =
    webbur::sandia_sgmga_size(d,&dim_weights[0],level,
			      &compute1DPoints[0], duplicate_tol, growth_rate,
                              &growth_rules[0]);

  Teuchos::Array<int> sparse_order(ntot*d);
  Teuchos::Array<int> sparse_index(ntot*d);
  Teuchos::Array<int> sparse_unique_index(num_total_pts);
  quad_points.resize(ntot);
  quad_weights.resize(ntot);
  quad_values.resize(ntot);
  Teuchos::Array<value_type> gp(ntot*d);

  webbur::sandia_sgmga_unique_index(d, &dim_weights[0], level,
				    &compute1DPoints[0],
				    duplicate_tol, ntot, num_total_pts,
				    growth_rate, &growth_rules[0],
				    &sparse_unique_index[0] );


  webbur::sandia_sgmga_index(d, &dim_weights[0], level,
			     ntot, num_total_pts, 
			     &sparse_unique_index[0],
			     growth_rate, &growth_rules[0],
			     &sparse_order[0], &sparse_index[0]);

  webbur::sandia_sgmga_weight(d,&dim_weights[0],level,
			      &compute1DWeights[0],
			      ntot, num_total_pts, &sparse_unique_index[0],
			      growth_rate, &growth_rules[0],
			      &quad_weights[0]);

  webbur::sandia_sgmga_point(d, &dim_weights[0], level,
			     &compute1DPoints[0],
			     ntot, &sparse_order[0], &sparse_index[0],
			     growth_rate, &growth_rules[0],
			     &gp[0]);

  for (ordinal_type i=0; i<ntot; i++) {
    quad_values[i].resize(sz);
    quad_points[i].resize(d);
    for (ordinal_type j=0; j<d; j++)
      quad_points[i][j] = gp[i*d+j];
    product_basis->evaluateBases(quad_points[i], quad_values[i]);
  }

  std::cout << "Number of quadrature points = " << ntot << std::endl;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::AnisoSparseGridQuadrature<ordinal_type, value_type>::
getQuadPoints() const
{
  return quad_points;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::AnisoSparseGridQuadrature<ordinal_type, value_type>::
getQuadWeights() const
{
  return quad_weights;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::AnisoSparseGridQuadrature<ordinal_type, value_type>::
getBasisAtQuadPoints() const
{
  return quad_values;
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::AnisoSparseGridQuadrature<ordinal_type,value_type>::
getMyPoints( int order, int dim, double x[] )
{
  Teuchos::Array<double> quad_points;
  Teuchos::Array<double> quad_weights;
  Teuchos::Array< Teuchos::Array<double> > quad_values;
  sgq->coordinate_bases[dim]->getQuadPoints(2*order-1, quad_points, 
                                            quad_weights, quad_values);
  for (int i = 0; i<quad_points.size(); i++) {
    x[i] = quad_points[i];
  }
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::AnisoSparseGridQuadrature<ordinal_type,value_type>::
getMyWeights( int order, int dim, double w[] )
{
  Teuchos::Array<double> quad_points;
  Teuchos::Array<double> quad_weights;
  Teuchos::Array< Teuchos::Array<double> > quad_values;
  sgq->coordinate_bases[dim]->getQuadPoints(2*order-1, quad_points, 
                                            quad_weights, quad_values);
  for (int i = 0; i<quad_points.size(); i++) {
    w[i] = quad_weights[i];
  }
}
