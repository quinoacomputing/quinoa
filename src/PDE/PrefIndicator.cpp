// *****************************************************************************
/*!
  \file      src/PDE/PrefIndicator.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Adaptive indicators for p-adaptive discontiunous Galerkin methods
  \details   This file contains functions that provide adaptive indicator
    function calculations for marking the number of degree of freedom of each
    element.
*/
// *****************************************************************************

#include "PrefIndicator.hpp"

#include "Tags.hpp"
#include "Vector.hpp"
#include "Integrate/Basis.hpp"
#include "Integrate/Quadrature.hpp"

namespace inciter {

static
void spectral_decay( std::size_t nunk,
                     const std::vector< int >& esuel,
                     const tk::Fields& unk,
                     std::size_t ndof,
                     std::size_t ndofmax,
                     tk::real tolref,
                     std::vector< std::size_t >& ndofel )
// *****************************************************************************
//! Evaluate the spectral-decay indicator and mark the ndof for each element
//! \param[in] nunk Number of unknowns
//! \param[in] esuel Elements surrounding elements
//! \param[in] unk Array of unknowns
//! \param[in] ndof Number of degrees of freedom in the solution
//! \param[in] ndofmax Max number of degrees of freedom for p-refinement
//! \param[in] tolref Tolerance for p-refinement
//! \param[in,out] ndofel Vector of local number of degrees of freedome
//! \details The spectral decay indicator, implemented in this functiopn,
//!   calculates the difference between the projections of the numerical
//!   solutions on finite element space of order p and p-1.
//! \see F. Naddei, et. al., "A comparison of refinement indicators for the
//!    p-adaptive simulation of steady and unsteady flows with discontinuous
//!    Galerkin methods" at https://doi.org/10.1016/j.jcp.2018.09.045, and G.
//!    Gassner, et al., "Explicit discontinuous Galerkin schemes with adaptation
//!    in space and time"
// *****************************************************************************
{
  const auto ncomp = unk.nprop() / ndof;

  // The array storing the adaptive indicator for each elements
  std::vector< tk::real > Ind(nunk, 0);

  for (std::size_t e=0; e<esuel.size()/4; ++e)
  {
    if(ndofel[e] > 1)
    {
      auto ng = tk::NGvol(ndofel[e]);

      // arrays for quadrature points
      std::array< std::vector< tk::real >, 3 > coordgp;
      std::vector< tk::real > wgp;

      coordgp[0].resize( ng );
      coordgp[1].resize( ng );
      coordgp[2].resize( ng );
      wgp.resize( ng );

      tk::GaussQuadratureTet( ng, coordgp, wgp );

      tk::real dU(0), U(0);

      // Gaussian quadrature
      for (std::size_t igp=0; igp<ng; ++igp)
      {
        // Compute the basis function
        auto B = tk::eval_basis( ndofel[e], coordgp[0][igp], coordgp[1][igp],
                                 coordgp[2][igp] );

        auto state = tk::eval_state( ncomp, 0, ndof, ndofel[e], e, unk, B );

        U += wgp[igp] * state[0] * state[0];

        if(ndofel[e] > 4)
        {
           auto dU_p2 = unk(e, 4, 0) * B[4]
                      + unk(e, 5, 0) * B[5]
                      + unk(e, 6, 0) * B[6]
                      + unk(e, 7, 0) * B[7]
                      + unk(e, 8, 0) * B[8]
                      + unk(e, 9, 0) * B[9];

           dU += wgp[igp] * dU_p2 * dU_p2;
        }
        else
        {
           auto dU_p1 = unk(e, 1, 0) * B[1]
                      + unk(e, 2, 0) * B[2]
                      + unk(e, 3, 0) * B[3];

           dU += wgp[igp] * dU_p1 * dU_p1;
        }
      }

      Ind[e] = log10( dU / U );
    }
  }

  // As for spectral-decay indicator, rho_p - rho_(p-1) actually is the leading
  // term of discretization error for the numerical solution of p-1. Therefore,
  // this function represents the qualitative behavior of the discretization
  // error. If the value is less than epsL which means the discretization error
  // is already a really small number, then the element should be de-refined. On
  // the other hand, if the value is larger than epsH which means the
  // discretization error is relatively large, then it should be refined.

  // Note: Spectral-decay indicator is a measurement of the continuity of the
  // numerical solution inside this element. So when this indicator appears
  // to be relatively large, there might be a shock inside this element and a
  // derefinement or h-refinement should be applied. This condition will be
  // implemented later.

  // As for the discretiazation-error based indicator, like spectral-decay
  // indicator, the choices for epsH and epsL are based on the data from
  // numerical experiments. Empirically, we found that when the epsH belongs
  // to [-4, -8] and epsL belongs to [-6, -14], decent results could be
  // observed. And then a linear projection is performed to map epsL and espH
  // to the range of [0, 1] so that it could be controlled by tolref.

  auto epsH = -4 - tolref * 4.0;
  auto epsL = -6 - tolref * 8.0;

  // Marke the ndof according to the adaptive indicator
  for (std::size_t e=0; e<esuel.size()/4; ++e)
  {
    if(Ind[e] < epsL)                          // Derefinement
    {
      if(ndofel[e] == 4)
        ndofel[e] = 1;
      else if(ndofel[e] == 10)
        ndofel[e] = 4;
    }
    else if(Ind[e] > epsH)                     // Refinement
    {
      if(ndofel[e] == 4 && ndofmax > 4)
        ndofel[e] = 10;
    }
  }
}

static
void non_conformity( std::size_t nunk,
                     std::size_t Nbfac,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     const std::vector< int >& esuel,
                     const std::vector< int >& esuf,
                     const std::vector< std::size_t >& inpofa,
                     const tk::Fields& unk,
                     std::size_t ndof,
                     std::size_t ndofmax,
                     std::vector< std::size_t >& ndofel )
// *****************************************************************************
//! Evaluate the non-conformity indicator and mark the ndof for each element
//! \param[in] nunk Number of unknowns
//! \param[in] Nbfac Number of internal faces
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] esuel Elements surrounding elements
//! \param[in] inpofa Face-node connectivity
//! \param[in] unk Array of unknowns
//! \param[in] ndof Number of degrees of freedom in the solution
//! \param[in] ndofmax Max number of degrees of freedom for p-refinement
//! \param[in,out] ndofel Vector of local number of degrees of freedome
//! \details The non-conformity indicator, this function implements, evaluates
//!   the jump in the numerical solutions as a measure of the numerical error.
//! \see F. Naddei, et. al., "A comparison of refinement indicators for the
//!   p-adaptive simulation of steady and unsteady flows with discontinuous
//!   Galerkin methods at https://doi.org/10.1016/j.jcp.2018.09.045.
//! \warning This indicator can only be applied in serial, i.e., single CPU, for
//!    now because the solution communication happens before eval_ndof() in DG,
//!    which will lead to incorrect evaluation of the numerical solution at the
//!    neighboring cells.
// *****************************************************************************
{
  const auto ncomp = unk.nprop() / ndof;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // The array storing the adaptive indicator for each elements
  std::vector< tk::real > Ind(nunk, 0);

  // compute error indicator for each face
  for (auto f=Nbfac; f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    auto ng_l = tk::NGfa(ndofel[el]);
    auto ng_r = tk::NGfa(ndofel[er]);

    // When the number of gauss points for the left and right element are
    // different, choose the larger ng
    auto ng = std::max( ng_l, ng_r );

    // arrays for quadrature points
    std::array< std::vector< tk::real >, 2 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    wgp.resize( ng );

    // get quadrature point weights and coordinates for triangle
    tk::GaussQuadratureTri( ng, coordgp, wgp );

    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel_l {{
      {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
      {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
      {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
      {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

    std::array< std::array< tk::real, 3>, 4 > coordel_r {{
      {{ cx[ inpoel[4*er  ] ], cy[ inpoel[4*er  ] ], cz[ inpoel[4*er  ] ] }},
      {{ cx[ inpoel[4*er+1] ], cy[ inpoel[4*er+1] ], cz[ inpoel[4*er+1] ] }},
      {{ cx[ inpoel[4*er+2] ], cy[ inpoel[4*er+2] ], cz[ inpoel[4*er+2] ] }},
      {{ cx[ inpoel[4*er+3] ], cy[ inpoel[4*er+3] ], cz[ inpoel[4*er+3] ] }} }};

    // Compute the determinant of Jacobian matrix
    auto detT_l =
      tk::Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    auto detT_r =
      tk::Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );

    // Extract the face coordinates
    std::array< std::array< tk::real, 3>, 3 > coordfa {{
      {{ cx[ inpofa[3*f  ] ], cy[ inpofa[3*f  ] ], cz[ inpofa[3*f  ] ] }},
      {{ cx[ inpofa[3*f+1] ], cy[ inpofa[3*f+1] ], cz[ inpofa[3*f+1] ] }},
      {{ cx[ inpofa[3*f+2] ], cy[ inpofa[3*f+2] ], cz[ inpofa[3*f+2] ] }} }};

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = tk::eval_gp( igp, coordfa, coordgp );

      //Compute the basis functions
      auto B_l = tk::eval_basis( ndofel[el],
        tk::Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
        tk::Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
        tk::Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l );
      auto B_r = tk::eval_basis( ndofel[er],
        tk::Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
        tk::Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
        tk::Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r );

      std::array< std::vector< tk::real >, 2 > state;

      state[0] = tk::eval_state( ncomp, 0, ndof, ndofel[el], el, unk, B_l );
      state[1] = tk::eval_state( ncomp, 0, ndof, ndofel[er], er, unk, B_r );

      Assert( unk[0].size() == ncomp, "Size mismatch" );
      Assert( unk[1].size() == ncomp, "Size mismatch" );

      auto rhoL = state[0][0];
      auto rhoR = state[1][0];

      auto ind = fabs( rhoL - rhoR ) / 2.0 * ( rhoL + rhoR );
      Ind[el] = std::max( ind, Ind[el] );
      Ind[er] = std::max( ind, Ind[er] );
    }
  }

  // By assuming a smooth solution, we use the non-conformity indicator to
  // represent the error for the numerical solution qualitatively. If the value
  // is less than epsL which means the error is already a really small number,
  // then the element should be de-refined. On the other hand, if the value is
  // larger than epsH which means the error is relatively large, then it should
  // be refined.

  // Marke the ndof according to the adaptive indicator
  for (std::size_t e=0; e<esuel.size()/4; ++e)
  {
    if(Ind[e] < 1e-4)                         // Derefinement
    {
      if(ndofel[e] == 10)
        ndofel[e] = 4;
      else if(ndofel[e] == 4)
        ndofel[e] = 1;
    }
    else if(Ind[e] > 1e-3)                    // Refinement
    {
      if(ndofel[e] == 4 && ndofmax > 4)
        ndofel[e] = 10;
      else if(ndofel[e] == 1)
        ndofel[e] = 4;
    }
  }
}

void eval_ndof( std::size_t nunk,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& inpoel,
                const inciter::FaceData& fd,
                const tk::Fields& unk,
                inciter::ctr::PrefIndicatorType indicator,
                std::size_t ndof,
                std::size_t ndofmax,
                tk::real tolref,
                std::vector< std::size_t >& ndofel )
// *****************************************************************************
//! Evaluate the adaptive indicator and mark the ndof for each element
//! \param[in] nunk Number of unknowns
//! \param[in] coord Array of nodal coordinates
//! \param[in] inpoel Element-node connectivity
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] unk Array of unknowns
//! \param[in] indicator p-refinement indicator type
//! \param[in] ndof Number of degrees of freedom in the solution
//! \param[in] ndofmax Max number of degrees of freedom for p-refinement
//! \param[in] tolref Tolerance for p-refinement
//! \param[in,out] ndofel Vector of local number of degrees of freedome
// *****************************************************************************
{
  const auto& esuel = fd.Esuel();

  if(indicator == inciter::ctr::PrefIndicatorType::SPECTRAL_DECAY)
    spectral_decay( nunk, esuel, unk, ndof, ndofmax, tolref, ndofel );
  else if(indicator == inciter::ctr::PrefIndicatorType::NON_CONFORMITY)
    non_conformity( nunk, fd.Nbfac(), inpoel, coord, esuel, fd.Esuf(),
                    fd.Inpofa(), unk, ndof, ndofmax, ndofel );
  else
    Throw( "No such adaptive indicator type" );
}

}
// inciter::
