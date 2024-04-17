// *****************************************************************************
/*!
  \file      src/PDE/Reconstruction.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Reconstruction for reconstructed discontinuous Galerkin methods
  \details   This file contains functions that reconstruct an "n"th order
    polynomial to an "n+1"th order polynomial using a least-squares
    reconstruction procedure.
*/
// *****************************************************************************

#include <array>
#include <vector>
#include <iostream>
#include <iomanip>

#include "Vector.hpp"
#include "Around.hpp"
#include "Base/HashMapReducer.hpp"
#include "Reconstruction.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Limiter.hpp"

namespace inciter {
extern ctr::InputDeck g_inputdeck;
}

namespace tk {

void
lhsLeastSq_P0P1( const inciter::FaceData& fd,
                 const Fields& geoElem,
                 const Fields& geoFace,
                 std::vector< std::array< std::array< real, 3 >, 3 > >& lhs_ls )
// *****************************************************************************
//  Compute lhs matrix for the least-squares reconstruction
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoElem Element geometry array
//! \param[in] geoFace Face geometry array
//! \param[in,out] lhs_ls LHS reconstruction matrix
//! \details This function computing the lhs matrix for reconstruction, is
//!   common for primitive and conserved quantities.
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();
  const auto nelem = fd.Esuel().size()/4;

  // Compute internal and boundary face contributions
  for (std::size_t f=0; f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1, "Left-side element detected as -1" );

    auto el = static_cast< std::size_t >(esuf[2*f]);
    auto er = esuf[2*f+1];

    std::array< real, 3 > geoElemR;
    std::size_t eR(0);

    // A second-order (piecewise linear) solution polynomial can be obtained
    // from the first-order (piecewise constant) FV solutions by using a
    // least-squares (LS) reconstruction process. LS uses the first-order
    // solutions from the cell being processed, and the cells surrounding it.
    // The LS system is obtaining by requiring the following to hold:
    // 'Taylor expansions of solution from cell-i to the centroids of each of
    // its neighboring cells should be equal to the cell average solution on
    // that neighbor cell.'
    // This gives a system of equations for the three second-order DOFs that are
    // to be determined. In 3D tetrahedral meshes, this would give four
    // equations (one for each neighbor )for the three unknown DOFs. This
    // overdetermined system is solved in the least-squares sense using the
    // normal equations approach. The normal equations approach involves
    // pre-multiplying the overdetermined system by the transpose of the system
    // matrix to obtain a square matrix (3x3 in this case).

    // get a 3x3 system by applying the normal equation approach to the
    // least-squares overdetermined system

    if (er > -1) {
    // internal face contribution
      eR = static_cast< std::size_t >(er);
      // Put in cell-centroid coordinates
      geoElemR = {{ geoElem(eR,1), geoElem(eR,2), geoElem(eR,3) }};
    }
    else {
    // boundary face contribution
      // Put in face-centroid coordinates
      geoElemR = {{ geoFace(f,4), geoFace(f,5), geoFace(f,6) }};
    }

    std::array< real, 3 > wdeltax{{ geoElemR[0]-geoElem(el,1),
                                    geoElemR[1]-geoElem(el,2),
                                    geoElemR[2]-geoElem(el,3) }};

    // define a lambda for contributing to lhs matrix
    auto lhs = [&]( std::size_t e ){
    for (std::size_t idir=0; idir<3; ++idir)
      for (std::size_t jdir=0; jdir<3; ++jdir)
        lhs_ls[e][idir][jdir] += wdeltax[idir] * wdeltax[jdir];
    };

    // always add left element contribution (at a boundary face, the internal
    // element is always the left element)
    lhs(el);
    // add right element contribution for internal faces only
    if (er > -1)
      if (eR < nelem) lhs(eR);

  }
}

void
intLeastSq_P0P1( const std::size_t rdof,
                 const inciter::FaceData& fd,
                 const Fields& geoElem,
                 const Fields& W,
                 std::vector< std::vector< std::array< real, 3 > > >& rhs_ls,
                 const std::vector< std::size_t >& varList )
// *****************************************************************************
//  \brief Compute internal surface contributions to rhs vector of the
//    least-squares reconstruction
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoElem Element geometry array
//! \param[in] W Solution vector to be reconstructed at recent time step
//! \param[in,out] rhs_ls RHS reconstruction vector
//! \param[in] varList List of indices in W, that need to be reconstructed
//! \details This function computing the internal face contributions to the rhs
//!   vector for reconstruction, is common for primitive and conserved
//!   quantities. If `W` == `U`, compute internal face contributions for the
//!   conserved variables. If `W` == `P`, compute internal face contributions
//!   for the primitive variables.
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();
  const auto nelem = fd.Esuel().size()/4;

  // Compute internal face contributions
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    auto el = static_cast< std::size_t >(esuf[2*f]);
    auto er = static_cast< std::size_t >(esuf[2*f+1]);

    // get a 3x3 system by applying the normal equation approach to the
    // least-squares overdetermined system

    // 'wdeltax' is the distance vector between the centroids of this element
    // and its neighbor
    std::array< real, 3 > wdeltax{{ geoElem(er,1)-geoElem(el,1),
                                    geoElem(er,2)-geoElem(el,2),
                                    geoElem(er,3)-geoElem(el,3) }};

    for (std::size_t idir=0; idir<3; ++idir)
    {
      // rhs vector
      for (std::size_t i=0; i<varList.size(); ++i)
      {
        auto c = varList[i];
        auto mark = c*rdof;
        rhs_ls[el][c][idir] +=
          wdeltax[idir] * (W(er,mark)-W(el,mark));
        if (er < nelem)
          rhs_ls[er][c][idir] +=
            wdeltax[idir] * (W(er,mark)-W(el,mark));
      }
    }
  }
}

void
bndLeastSqConservedVar_P0P1(
  ncomp_t ncomp,
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  real t,
  const StateFn& state,
  const Fields& P,
  const Fields& U,
  std::vector< std::vector< std::array< real, 3 > > >& rhs_ls,
  const std::vector< std::size_t >& varList,
  std::size_t nprim )
// *****************************************************************************
//  \brief Compute boundary surface contributions to rhs vector of the
//    least-squares reconstruction of conserved quantities of the PDE system
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] t Physical time
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] P Primitive vector to be reconstructed at recent time step
//! \param[in] U Solution vector to be reconstructed at recent time step
//! \param[in,out] rhs_ls RHS reconstruction vector
//! \param[in] varList List of indices in W, that need to be reconstructed
//! \param[in] nprim This is the number of primitive quantities stored for this
//!   PDE system. This is necessary to extend the state vector to the right
//!   size, so that correct boundary conditions are obtained.
//!   A default is set to 0, so that calling code for systems that do not store
//!   primitive quantities does not need to specify this argument.
//! \details This function computing the boundary face contributions to the rhs
//!   vector for reconstruction, is used for conserved quantities only.
// *****************************************************************************
{
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find(static_cast<int>(s));// faces for side set
    if (bc != end(bface))
    {
      // Compute boundary face contributions
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "physical boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);

        // arrays for quadrature points
        std::array< real, 3 >
          fc{{ geoFace(f,4), geoFace(f,5), geoFace(f,6) }};
        std::array< real, 3 >
          fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

        // Compute the state variables at the left element
        std::vector< real >B(1,1.0);
        auto ul = eval_state( ncomp, rdof, 1, el, U, B );
        auto uprim = eval_state( nprim, rdof, 1, el, P, B );

        // consolidate primitives into state vector
        ul.insert(ul.end(), uprim.begin(), uprim.end());

        Assert( ul.size() == ncomp+nprim, "Incorrect size for "
                "appended state vector" );

        // Compute the state at the face-center using BC
        auto ustate = state( ncomp, mat_blk, ul, fc[0], fc[1], fc[2], t, fn );

        std::array< real, 3 > wdeltax{{ fc[0]-geoElem(el,1),
                                        fc[1]-geoElem(el,2),
                                        fc[2]-geoElem(el,3) }};

        for (std::size_t idir=0; idir<3; ++idir)
        {
          // rhs vector
          for (std::size_t i=0; i<varList.size(); ++i)
          {
            auto c = varList[i];
            rhs_ls[el][c][idir] +=
              wdeltax[idir] * (ustate[1][c]-ustate[0][c]);
          }
        }
      }
    }
  }
}

void
solveLeastSq_P0P1(
  const std::size_t rdof,
  const std::vector< std::array< std::array< real, 3 >, 3 > >& lhs,
  const std::vector< std::vector< std::array< real, 3 > > >& rhs,
  Fields& W,
  const std::vector< std::size_t >& varList )
// *****************************************************************************
//  Solve the 3x3 linear system for least-squares reconstruction
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] lhs LHS reconstruction matrix
//! \param[in] rhs RHS reconstruction vector
//! \param[in,out] W Solution vector to be reconstructed at recent time step
//! \param[in] varList List of indices in W, that need to be reconstructed
//! \details Solves the 3x3 linear system for each element, individually. For
//!   systems that require reconstructions of primitive quantities, this should
//!   be called twice, once with the argument 'W' as U (conserved), and again
//!   with 'W' as P (primitive).
// *****************************************************************************
{
  auto nelem = lhs.size();

  for (std::size_t e=0; e<nelem; ++e)
  {
    for (std::size_t i=0; i<varList.size(); ++i)
    {
      auto mark = varList[i]*rdof;

      // solve system using Cramer's rule
      auto ux = tk::cramer( lhs[e], rhs[e][varList[i]] );

      W(e,mark+1) = ux[0];
      W(e,mark+2) = ux[1];
      W(e,mark+3) = ux[2];
    }
  }
}

void
recoLeastSqExtStencil(
  std::size_t rdof,
  std::size_t e,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const Fields& geoElem,
  Fields& W,
  const std::vector< std::size_t >& varList )
// *****************************************************************************
//  \brief Reconstruct the second-order solution using least-squares approach
//    from an extended stencil involving the node-neighbors
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] e Element whoes solution is being reconstructed
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element-node connectivity
//! \param[in] geoElem Element geometry array
//! \param[in,out] W Solution vector to be reconstructed at recent time step
//! \param[in] varList List of indices in W, that need to be reconstructed
//! \details A second-order (piecewise linear) solution polynomial is obtained
//!   from the first-order (piecewise constant) FV solutions by using a
//!   least-squares (LS) reconstruction process. This LS reconstruction function
//!   using the nodal-neighbors of a cell, to get an overdetermined system of
//!   equations for the derivatives of the solution. This overdetermined system
//!   is solved in the least-squares sense using the normal equations approach.
// *****************************************************************************
{
  // lhs matrix
  std::array< std::array< tk::real, 3 >, 3 >
    lhs_ls( {{ {{0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0}} }} );
  // rhs matrix
  std::vector< std::array< tk::real, 3 > >
  rhs_ls( varList.size(), {{ 0.0, 0.0, 0.0 }} );

  // loop over all nodes of the element e
  for (std::size_t lp=0; lp<4; ++lp)
  {
    auto p = inpoel[4*e+lp];
    const auto& pesup = cref_find(esup, p);

    // loop over all the elements surrounding this node p
    for (auto er : pesup)
    {
      // centroid distance
      std::array< real, 3 > wdeltax{{ geoElem(er,1)-geoElem(e,1),
                                      geoElem(er,2)-geoElem(e,2),
                                      geoElem(er,3)-geoElem(e,3) }};

      // contribute to lhs matrix
      for (std::size_t idir=0; idir<3; ++idir)
        for (std::size_t jdir=0; jdir<3; ++jdir)
          lhs_ls[idir][jdir] += wdeltax[idir] * wdeltax[jdir];

      // compute rhs matrix
      for (std::size_t i=0; i<varList.size(); i++)
      {
        auto mark = varList[i]*rdof;
        for (std::size_t idir=0; idir<3; ++idir)
          rhs_ls[i][idir] +=
            wdeltax[idir] * (W(er,mark)-W(e,mark));

      }
    }
  }

  // solve least-square normal equation system using Cramer's rule
  for (std::size_t i=0; i<varList.size(); i++)
  {
    auto mark = varList[i]*rdof;

    auto ux = tk::cramer( lhs_ls, rhs_ls[i] );

    // Update the P1 dofs with the reconstructioned gradients.
    // Since this reconstruction does not affect the cell-averaged solution,
    // W(e,mark+0) is unchanged.
    W(e,mark+1) = ux[0];
    W(e,mark+2) = ux[1];
    W(e,mark+3) = ux[2];
  }
}

void
transform_P0P1( std::size_t rdof,
                std::size_t nelem,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                Fields& W,
                const std::vector< std::size_t >& varList )
// *****************************************************************************
//  Transform the reconstructed P1-derivatives to the Dubiner dofs
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] W Second-order reconstructed vector which gets transformed to
//!   the Dubiner reference space
//! \param[in] varList List of indices in W, that need to be reconstructed
//! \details Since the DG solution (and the primitive quantities) are assumed to
//!   be stored in the Dubiner space, this transformation from Taylor
//!   coefficients to Dubiner coefficients is necessary.
// *****************************************************************************
{
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<nelem; ++e)
  {
    // Extract the element coordinates
    std::array< std::array< real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
    }};

    auto jacInv =
      tk::inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    // Compute the derivatives of basis function for DG(P1)
    auto dBdx = tk::eval_dBdx_p1( rdof, jacInv );

    for (std::size_t i=0; i<varList.size(); ++i)
    {
      auto mark = varList[i]*rdof;

      // solve system using Cramer's rule
      auto ux = tk::cramer( {{ {{dBdx[0][1], dBdx[0][2], dBdx[0][3]}},
                               {{dBdx[1][1], dBdx[1][2], dBdx[1][3]}},
                               {{dBdx[2][1], dBdx[2][2], dBdx[2][3]}} }},
                            {{ W(e,mark+1),
                               W(e,mark+2),
                               W(e,mark+3) }} );

      // replace physical derivatives with transformed dofs
      W(e,mark+1) = ux[0];
      W(e,mark+2) = ux[1];
      W(e,mark+3) = ux[2];
    }
  }
}

void
THINCReco( std::size_t rdof,
           std::size_t nmat,
           std::size_t e,
           const std::vector< std::size_t >& inpoel,
           const UnsMesh::Coords& coord,
           const Fields& geoElem,
           const std::array< real, 3 >& ref_xp,
           const Fields& U,
           const Fields& P,
           bool intInd,
           const std::vector< std::size_t >& matInt,
           [[maybe_unused]] const std::vector< real >& vfmin,
           [[maybe_unused]] const std::vector< real >& vfmax,
           std::vector< real >& state )
// *****************************************************************************
//  Compute THINC reconstructions at quadrature point for multi-material flows
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nmat Total number of materials
//! \param[in] e Element for which interface reconstruction is being calculated
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] ref_xp Quadrature point in reference space
//! \param[in] U Solution vector
//! \param[in] P Vector of primitives
//! \param[in] intInd Boolean which indicates if the element contains a
//!   material interface
//! \param[in] matInt Array indicating which material has an interface
//! \param[in] vfmin Vector containing min volume fractions for each material
//!   in this cell
//! \param[in] vfmax Vector containing max volume fractions for each material
//!   in this cell
//! \param[in,out] state Unknown/state vector at quadrature point, modified
//!   if near interfaces using THINC
//! \details This function is an interface for the multimat PDEs that use the
//!   algebraic multi-material THINC reconstruction. This particular function
//!   should only be called for multimat.
// *****************************************************************************
{
  using inciter::volfracDofIdx;
  using inciter::densityDofIdx;
  using inciter::momentumDofIdx;
  using inciter::energyDofIdx;
  using inciter::pressureDofIdx;
  using inciter::velocityDofIdx;
  using inciter::deformDofIdx;
  using inciter::stressDofIdx;
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::pressureIdx;
  using inciter::velocityIdx;
  using inciter::deformIdx;
  using inciter::stressIdx;

  auto bparam = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp_param >();
  const auto ncomp = U.nprop()/rdof;
  const auto& solidx = inciter::g_inputdeck.get< tag::matidxmap,
    tag::solidx >();

  // Step-1: Perform THINC reconstruction
  // create a vector of volume-fractions and pass it to the THINC function
  std::vector< real > alSol(rdof*nmat, 0.0);
  std::vector< real > alReco(nmat, 0.0);
  for (std::size_t k=0; k<nmat; ++k) {
    auto mark = k*rdof;
    for (std::size_t i=0; i<rdof; ++i) {
      alSol[mark+i] = U(e, volfracDofIdx(nmat,k,rdof,i));
    }
    // initialize with TVD reconstructions which will be modified if near
    // material interface
    alReco[k] = state[volfracIdx(nmat,k)];
  }
  THINCFunction(rdof, nmat, e, inpoel, coord, ref_xp, geoElem(e,0), bparam,
    alSol, intInd, matInt, alReco);

  // check reconstructed volfracs for positivity
  bool neg_vf = false;
  for (std::size_t k=0; k<nmat; ++k) {
    if (alReco[k] < 1e-16 && matInt[k] > 0) neg_vf = true;
  }
  for (std::size_t k=0; k<nmat; ++k) {
    if (neg_vf) {
      std::cout << "Material-id:        " << k << std::endl;
      std::cout << "Volume-fraction:    " << std::setprecision(18) << alReco[k]
        << std::endl;
      std::cout << "Cell-avg vol-frac:  " << std::setprecision(18) <<
        U(e,volfracDofIdx(nmat,k,rdof,0)) << std::endl;
      std::cout << "Material-interface? " << intInd << std::endl;
      std::cout << "Mat-k-involved?     " << matInt[k] << std::endl;
    }
  }
  if (neg_vf) Throw("Material has negative volume fraction after THINC "
    "reconstruction.");

  // Step-2: Perform consistent reconstruction on other conserved quantities
  if (intInd)
  {
    auto rhobCC(0.0), rhobHO(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alCC = U(e, volfracDofIdx(nmat,k,rdof,0));
      alCC = std::max(1e-14, alCC);

      if (matInt[k])
      {
        state[volfracIdx(nmat,k)] = alReco[k];
        state[densityIdx(nmat,k)] = alReco[k]
          * U(e, densityDofIdx(nmat,k,rdof,0))/alCC;
        state[energyIdx(nmat,k)] = alReco[k]
          * U(e, energyDofIdx(nmat,k,rdof,0))/alCC;
        state[ncomp+pressureIdx(nmat,k)] = alReco[k]
          * P(e, pressureDofIdx(nmat,k,rdof,0))/alCC;
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              state[deformIdx(nmat,solidx[k],i,j)] =
                U(e, deformDofIdx(nmat,solidx[k],i,j,rdof,0));

          for (std::size_t i=0; i<6; ++i)
            state[ncomp+stressIdx(nmat,solidx[k],i)] = alReco[k]
              * P(e, stressDofIdx(nmat,solidx[k],i,rdof,0))/alCC;
        }
      }

      rhobCC += U(e, densityDofIdx(nmat,k,rdof,0));
      rhobHO += state[densityIdx(nmat,k)];
    }

    // consistent reconstruction for bulk momentum
    for (std::size_t i=0; i<3; ++i)
    {
      state[momentumIdx(nmat,i)] = rhobHO
        * U(e, momentumDofIdx(nmat,i,rdof,0))/rhobCC;
      state[ncomp+velocityIdx(nmat,i)] =
        P(e, velocityDofIdx(nmat,i,rdof,0));
    }
  }
}

void
THINCRecoTransport( std::size_t rdof,
                    std::size_t,
                    std::size_t e,
                    const std::vector< std::size_t >& inpoel,
                    const UnsMesh::Coords& coord,
                    const Fields& geoElem,
                    const std::array< real, 3 >& ref_xp,
                    const Fields& U,
                    const Fields&,
                    [[maybe_unused]] const std::vector< real >& vfmin,
                    [[maybe_unused]] const std::vector< real >& vfmax,
                    std::vector< real >& state )
// *****************************************************************************
//  Compute THINC reconstructions at quadrature point for transport
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] e Element for which interface reconstruction is being calculated
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] ref_xp Quadrature point in reference space
//! \param[in] U Solution vector
//! \param[in] vfmin Vector containing min volume fractions for each material
//!   in this cell
//! \param[in] vfmax Vector containing max volume fractions for each material
//!   in this cell
//! \param[in,out] state Unknown/state vector at quadrature point, modified
//!   if near interfaces using THINC
//! \details This function is an interface for the transport PDEs that use the
//!   algebraic multi-material THINC reconstruction. This particular function
//!   should only be called for transport.
// *****************************************************************************
{
  auto bparam = inciter::g_inputdeck.get< tag::transport,
    tag::intsharp_param >();
  auto ncomp = U.nprop()/rdof;

  // interface detection
  std::vector< std::size_t > matInt(ncomp, 0);
  std::vector< tk::real > alAvg(ncomp, 0.0);
  for (std::size_t k=0; k<ncomp; ++k)
    alAvg[k] = U(e, k*rdof);
  auto intInd = inciter::interfaceIndicator(ncomp, alAvg, matInt);

  // create a vector of volume-fractions and pass it to the THINC function
  std::vector< real > alSol(rdof*ncomp, 0.0);
  // initialize with TVD reconstructions (modified if near interface)
  auto alReco = state;
  for (std::size_t k=0; k<ncomp; ++k) {
    auto mark = k*rdof;
    for (std::size_t i=0; i<rdof; ++i) {
      alSol[mark+i] = U(e,mark+i);
    }
  }
  THINCFunction(rdof, ncomp, e, inpoel, coord, ref_xp, geoElem(e,0), bparam,
    alSol, intInd, matInt, alReco);

  state = alReco;
}

void
THINCFunction( std::size_t rdof,
               std::size_t nmat,
               std::size_t e,
               const std::vector< std::size_t >& inpoel,
               const UnsMesh::Coords& coord,
               const std::array< real, 3 >& ref_xp,
               real vol,
               real bparam,
               const std::vector< real >& alSol,
               bool intInd,
               const std::vector< std::size_t >& matInt,
               std::vector< real >& alReco )
// *****************************************************************************
//  Old version of the Multi-Medium THINC reconstruction function
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nmat Total number of materials
//! \param[in] e Element for which interface reconstruction is being calculated
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] ref_xp Quadrature point in reference space
//! \param[in] vol Element volume
//! \param[in] bparam User specified Beta for THINC, from the input file
//! \param[in] alSol Volume fraction solution vector for element e
//! \param[in] intInd Interface indicator, true if e is interface element
//! \param[in] matInt Vector indicating materials which constitute interface
//! \param[in,out] alReco Unknown/state vector at quadrature point, which gets
//!   modified if near interface using MM-THINC
//! \details This function computes the interface reconstruction using the
//!   algebraic multi-material THINC reconstruction for each material at the
//!   given (ref_xp) quadrature point. This function is based on the following:
//!   Pandare A. K., Waltz J., & Bakosi J. (2021) Multi-Material Hydrodynamics
//!   with Algebraic Sharp Interface Capturing. Computers & Fluids,
//!   doi: https://doi.org/10.1016/j.compfluid.2020.104804.
//!   This function will be removed after the newer version (see
//!   THINCFunction_new) is sufficiently tested.
// *****************************************************************************
{
  // determine number of materials with interfaces in this cell
  auto epsl(1e-4), epsh(1e-1), bred(1.25), bmod(bparam);
  std::size_t nIntMat(0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    auto alk = alSol[k*rdof];
    if (alk > epsl)
    {
      ++nIntMat;
      if ((alk > epsl) && (alk < epsh))
        bmod = std::min(bmod,
          (alk-epsl)/(epsh-epsl) * (bred - bparam) + bparam);
      else if (alk > epsh)
        bmod = bred;
    }
  }

  if (nIntMat > 2) bparam = bmod;

  // compression parameter
  auto beta = bparam/std::cbrt(6.0*vol);

  if (intInd)
  {
    // 1. Get unit normals to material interface

    // Compute Jacobian matrix for converting Dubiner dofs to derivatives
    const auto& cx = coord[0];
    const auto& cy = coord[1];
    const auto& cz = coord[2];

    std::array< std::array< real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
    }};

    auto jacInv =
      tk::inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    auto dBdx = tk::eval_dBdx_p1( rdof, jacInv );

    std::array< real, 3 > nInt;
    std::vector< std::array< real, 3 > > ref_n(nmat, {{0.0, 0.0, 0.0}});

    // Get normals
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Get derivatives from moments in Dubiner space
      for (std::size_t i=0; i<3; ++i)
        nInt[i] = dBdx[i][1] * alSol[k*rdof+1]
          + dBdx[i][2] * alSol[k*rdof+2]
          + dBdx[i][3] * alSol[k*rdof+3];

      auto nMag = std::sqrt(tk::dot(nInt, nInt)) + 1e-14;

      for (std::size_t i=0; i<3; ++i)
        nInt[i] /= nMag;

      // project interface normal onto local/reference coordinate system
      for (std::size_t i=0; i<3; ++i)
      {
        std::array< real, 3 > axis{
          coordel[i+1][0]-coordel[0][0],
          coordel[i+1][1]-coordel[0][1],
          coordel[i+1][2]-coordel[0][2] };
        ref_n[k][i] = tk::dot(nInt, axis);
      }
    }

    // 2. Reconstruct volume fractions using THINC
    auto max_lim = 1.0 - (static_cast<tk::real>(nmat-1)*1.0e-12);
    auto min_lim = 1e-12;
    auto sum_inter(0.0), sum_non_inter(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      if (matInt[k])
      {
        // get location of material interface (volume fraction 0.5) from the
        // assumed tanh volume fraction distribution, and cell-averaged
        // volume fraction
        auto alCC(alSol[k*rdof]);
        auto Ac(0.0), Bc(0.0), Qc(0.0);
        if ((std::abs(ref_n[k][0]) > std::abs(ref_n[k][1]))
          && (std::abs(ref_n[k][0]) > std::abs(ref_n[k][2])))
        {
          Ac = std::exp(0.5*beta*ref_n[k][0]);
          Bc = std::exp(0.5*beta*(ref_n[k][1]+ref_n[k][2]));
          Qc = std::exp(0.5*beta*ref_n[k][0]*(2.0*alCC-1.0));
        }
        else if ((std::abs(ref_n[k][1]) > std::abs(ref_n[k][0]))
          && (std::abs(ref_n[k][1]) > std::abs(ref_n[k][2])))
        {
          Ac = std::exp(0.5*beta*ref_n[k][1]);
          Bc = std::exp(0.5*beta*(ref_n[k][0]+ref_n[k][2]));
          Qc = std::exp(0.5*beta*ref_n[k][1]*(2.0*alCC-1.0));
        }
        else
        {
          Ac = std::exp(0.5*beta*ref_n[k][2]);
          Bc = std::exp(0.5*beta*(ref_n[k][0]+ref_n[k][1]));
          Qc = std::exp(0.5*beta*ref_n[k][2]*(2.0*alCC-1.0));
        }
        auto d = std::log((1.0-Ac*Qc) / (Ac*Bc*(Qc-Ac))) / (2.0*beta);

        // THINC reconstruction
        auto al_c = 0.5 * (1.0 + std::tanh(beta*(tk::dot(ref_n[k], ref_xp) + d)));

        alReco[k] = std::min(max_lim, std::max(min_lim, al_c));

        sum_inter += alReco[k];
      } else
      {
        sum_non_inter += alReco[k];
      }
      // else, if this material does not have an interface close-by, the TVD
      // reconstructions must be used for state variables. This is ensured by
      // initializing the alReco vector as the TVD state.
    }

    // Rescale volume fractions of interface-materials to ensure unit sum
    auto sum_rest = 1.0 - sum_non_inter;
    for (std::size_t k=0; k<nmat; ++k)
      if(matInt[k])
        alReco[k] = alReco[k] * sum_rest / sum_inter;
  }
}

void
THINCFunction_new( std::size_t rdof,
                   std::size_t nmat,
                   std::size_t e,
                   const std::vector< std::size_t >& inpoel,
                   const UnsMesh::Coords& coord,
                   const std::array< real, 3 >& ref_xp,
                   real vol,
                   real bparam,
                   const std::vector< real >& alSol,
                   bool intInd,
                   const std::vector< std::size_t >& matInt,
                   std::vector< real >& alReco )
// *****************************************************************************
//  New Multi-Medium THINC reconstruction function for volume fractions
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nmat Total number of materials
//! \param[in] e Element for which interface reconstruction is being calculated
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] ref_xp Quadrature point in reference space
//! \param[in] vol Element volume
//! \param[in] bparam User specified Beta for THINC, from the input file
//! \param[in] alSol Volume fraction solution vector for element e
//! \param[in] intInd Interface indicator, true if e is interface element
//! \param[in] matInt Vector indicating materials which constitute interface
//! \param[in,out] alReco Unknown/state vector at quadrature point, which gets
//!   modified if near interface using MM-THINC
//! \details This function computes the interface reconstruction using the
//!   algebraic multi-material THINC reconstruction for each material at the
//!   given (ref_xp) quadrature point. This function succeeds the older version
//!   of the mm-THINC (see THINCFunction), but is still under testing and is
//!   currently experimental.
// *****************************************************************************
{
  // compression parameter
  auto beta = bparam/std::cbrt(6.0*vol);

  // If the cell is not material interface, return this function
  if (not intInd) return;

  // If the cell is material interface, THINC reconstruction is applied
  // Step 1. Get unit normals to material interface
  // -------------------------------------------------------------------------

  // Compute Jacobian matrix for converting Dubiner dofs to derivatives
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  std::array< std::array< real, 3>, 4 > coordel {{
    {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
    {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
    {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
    {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
  }};

  auto jacInv =
    tk::inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

  auto dBdx = tk::eval_dBdx_p1( rdof, jacInv );

  std::array< real, 3 > nInt;
  std::array< real, 3 > ref_n{0.0, 0.0, 0.0};
  auto almax(0.0);
  std::size_t kmax(0);

  // Determine index of material present in majority
  for (std::size_t k=0; k<nmat; ++k)
  {
    auto alk = alSol[k*rdof];
    if (alk > almax)
    {
      almax = alk;
      kmax = k;
    }
  }

  // Get normals of material present in majority
  // Get derivatives from moments in Dubiner space
  for (std::size_t i=0; i<3; ++i)
    nInt[i] = dBdx[i][1] * alSol[kmax*rdof+1]
      + dBdx[i][2] * alSol[kmax*rdof+2]
      + dBdx[i][3] * alSol[kmax*rdof+3];

  auto nMag = std::sqrt(tk::dot(nInt, nInt)) + 1e-14;

  for (std::size_t i=0; i<3; ++i)
    nInt[i] /= nMag;

  // project interface normal onto local/reference coordinate system
  for (std::size_t i=0; i<3; ++i)
  {
    std::array< real, 3 > axis{
      coordel[i+1][0]-coordel[0][0],
      coordel[i+1][1]-coordel[0][1],
      coordel[i+1][2]-coordel[0][2] };
    ref_n[i] = tk::dot(nInt, axis);
  }

  // Step 2. Reconstruct volume fraction of majority material using THINC
  // -------------------------------------------------------------------------

  auto al_max = 1.0 - (static_cast<tk::real>(nmat-1)*1.0e-12);
  auto al_min = 1e-12;
  auto alsum(0.0);
  // get location of material interface (volume fraction 0.5) from the
  // assumed tanh volume fraction distribution, and cell-averaged
  // volume fraction
  auto alCC(alSol[kmax*rdof]);
  auto Ac(0.0), Bc(0.0), Qc(0.0);
  if ((std::abs(ref_n[0]) > std::abs(ref_n[1]))
    && (std::abs(ref_n[0]) > std::abs(ref_n[2])))
  {
    Ac = std::exp(0.5*beta*ref_n[0]);
    Bc = std::exp(0.5*beta*(ref_n[1]+ref_n[2]));
    Qc = std::exp(0.5*beta*ref_n[0]*(2.0*alCC-1.0));
  }
  else if ((std::abs(ref_n[1]) > std::abs(ref_n[0]))
    && (std::abs(ref_n[1]) > std::abs(ref_n[2])))
  {
    Ac = std::exp(0.5*beta*ref_n[1]);
    Bc = std::exp(0.5*beta*(ref_n[0]+ref_n[2]));
    Qc = std::exp(0.5*beta*ref_n[1]*(2.0*alCC-1.0));
  }
  else
  {
    Ac = std::exp(0.5*beta*ref_n[2]);
    Bc = std::exp(0.5*beta*(ref_n[0]+ref_n[1]));
    Qc = std::exp(0.5*beta*ref_n[2]*(2.0*alCC-1.0));
  }
  auto d = std::log((1.0-Ac*Qc) / (Ac*Bc*(Qc-Ac))) / (2.0*beta);

  // THINC reconstruction
  auto al_c = 0.5 * (1.0 + std::tanh(beta*(tk::dot(ref_n, ref_xp) + d)));

  alReco[kmax] = std::min(al_max, std::max(al_min, al_c));
  alsum += alReco[kmax];

  // if this material does not have an interface close-by, the TVD
  // reconstructions must be used for state variables. This is ensured by
  // initializing the alReco vector as the TVD state.
  for (std::size_t k=0; k<nmat; ++k) {
    if (!matInt[k]) {
      alsum += alReco[k];
    }
  }

  // Step 3. Do multimaterial cell corrections
  // -------------------------------------------------------------------------

  // distribute remaining volume to rest of materials
  auto sum_left = 1.0 - alsum;
  real den = 0.0;
  for (std::size_t k=0; k<nmat; ++k) {
    if (matInt[k] && k != kmax) {
      auto mark = k * rdof;
      alReco[k] = sum_left * alSol[mark];
      den += alSol[mark];
    }
  }
  // the distributed volfracs might be below al_min, correct that
  real err = 0.0;
  for (std::size_t k=0; k<nmat; ++k) {
    if (matInt[k] && k != kmax) {
      alReco[k] /= den;
      if (alReco[k] < al_min) {
        err += al_min - alReco[k];
        alReco[k] = al_min;
      }
    }
  }

  // balance out errors
  alReco[kmax] -= err;
}

std::vector< tk::real >
evalPolynomialSol( const std::vector< inciter::EOS >& mat_blk,
                   int intsharp,
                   std::size_t ncomp,
                   std::size_t nprim,
                   std::size_t rdof,
                   std::size_t nmat,
                   std::size_t e,
                   std::size_t dof_e,
                   const std::vector< std::size_t >& inpoel,
                   const UnsMesh::Coords& coord,
                   const Fields& geoElem,
                   const std::array< real, 3 >& ref_gp,
                   const std::vector< real >& B,
                   const Fields& U,
                   const Fields& P )
// *****************************************************************************
//  Evaluate polynomial solution at quadrature point
//! \param[in] mat_blk EOS material block
//! \param[in] intsharp Interface reconstruction indicator
//! \param[in] ncomp Number of components in the PDE system
//! \param[in] nprim Number of primitive quantities
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nmat Total number of materials
//! \param[in] e Element for which polynomial solution is being evaluated
//! \param[in] dof_e Degrees of freedom for element
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] ref_gp Quadrature point in reference space
//! \param[in] B Basis function at given quadrature point
//! \param[in] U Solution vector
//! \param[in] P Vector of primitives
//! \return High-order unknown/state vector at quadrature point, modified
//!   if near interfaces using THINC
// *****************************************************************************
{
  std::vector< real > state;
  std::vector< real > sprim;

  state = eval_state( ncomp, rdof, dof_e, e, U, B );
  sprim = eval_state( nprim, rdof, dof_e, e, P, B );

  // interface detection
  std::vector< std::size_t > matInt(nmat, 0);
  bool intInd(false);
  if (nmat > 1) {
    std::vector< tk::real > alAvg(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
      alAvg[k] = U(e, inciter::volfracDofIdx(nmat,k,rdof,0));
    intInd = inciter::interfaceIndicator(nmat, alAvg, matInt);
  }

  // consolidate primitives into state vector
  state.insert(state.end(), sprim.begin(), sprim.end());

  if (intsharp > 0)
  {
    std::vector< tk::real > vfmax(nmat, 0.0), vfmin(nmat, 0.0);

    // Until the appropriate setup for activating THINC with Transport
    // is ready, the following two chunks of code will need to be commented
    // for using THINC with Transport
    //for (std::size_t k=0; k<nmat; ++k) {
    //  vfmin[k] = VolFracMax(el, 2*k, 0);
    //  vfmax[k] = VolFracMax(el, 2*k+1, 0);
    //}
    tk::THINCReco(rdof, nmat, e, inpoel, coord, geoElem,
      ref_gp, U, P, intInd, matInt, vfmin, vfmax, state);

    // Until the appropriate setup for activating THINC with Transport
    // is ready, the following lines will need to be uncommented for
    // using THINC with Transport
    //tk::THINCRecoTransport(rdof, nmat, el, inpoel, coord,
    //  geoElem, ref_gp_l, U, P, vfmin, vfmax, state[0]);
  }

  // physical constraints
  if (state.size() > ncomp) {
    using inciter::pressureIdx;
    using inciter::volfracIdx;
    using inciter::densityIdx;

    for (std::size_t k=0; k<nmat; ++k) {
      state[ncomp+pressureIdx(nmat,k)] = constrain_pressure( mat_blk,
        state[ncomp+pressureIdx(nmat,k)], state[densityIdx(nmat,k)],
        state[volfracIdx(nmat,k)], k );
    }
  }

  return state;
}

std::vector< tk::real >
evalFVSol( const std::vector< inciter::EOS >& mat_blk,
           int intsharp,
           std::size_t ncomp,
           std::size_t nprim,
           std::size_t rdof,
           std::size_t nmat,
           std::size_t e,
           const std::vector< std::size_t >& inpoel,
           const UnsMesh::Coords& coord,
           const Fields& geoElem,
           const std::array< real, 3 >& ref_gp,
           const std::vector< real >& B,
           const Fields& U,
           const Fields& P,
           int srcFlag )
// *****************************************************************************
//  Evaluate second-order FV solution at quadrature point
//! \param[in] mat_blk EOS material block
//! \param[in] intsharp Interface reconstruction indicator
//! \param[in] ncomp Number of components in the PDE system
//! \param[in] nprim Number of primitive quantities
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nmat Total number of materials
//! \param[in] e Element for which polynomial solution is being evaluated
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] ref_gp Quadrature point in reference space
//! \param[in] B Basis function at given quadrature point
//! \param[in] U Solution vector
//! \param[in] P Vector of primitives
//! \param[in] srcFlag Whether the energy source was added to element e
//! \return High-order unknown/state vector at quadrature point, modified
//!   if near interfaces using THINC
// *****************************************************************************
{
  using inciter::pressureIdx;
  using inciter::velocityIdx;
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::energyIdx;
  using inciter::momentumIdx;

  std::vector< real > state;
  std::vector< real > sprim;

  state = eval_state( ncomp, rdof, rdof, e, U, B );
  sprim = eval_state( nprim, rdof, rdof, e, P, B );

  // interface detection so that eos is called on the appropriate quantities
  std::vector< std::size_t > matInt(nmat, 0);
  std::vector< tk::real > alAvg(nmat, 0.0);
  for (std::size_t k=0; k<nmat; ++k)
    alAvg[k] = U(e, inciter::volfracDofIdx(nmat,k,rdof,0));
  auto intInd = inciter::interfaceIndicator(nmat, alAvg, matInt);

  // get mat-energy from reconstructed mat-pressure
  auto rhob(0.0);
  for (std::size_t k=0; k<nmat; ++k) {
    auto alk = state[volfracIdx(nmat,k)];
    if (matInt[k]) {
      alk = std::max(std::min(alk, 1.0-static_cast<tk::real>(nmat-1)*1e-12),
        1e-12);
    }
    state[energyIdx(nmat,k)] = alk *
      mat_blk[k].compute< inciter::EOS::totalenergy >(
      state[densityIdx(nmat,k)]/alk, sprim[velocityIdx(nmat,0)],
      sprim[velocityIdx(nmat,1)], sprim[velocityIdx(nmat,2)],
      sprim[pressureIdx(nmat,k)]/alk);
    rhob += state[densityIdx(nmat,k)];
  }
  // get momentum from reconstructed velocity and bulk density
  for (std::size_t i=0; i<3; ++i) {
    state[momentumIdx(nmat,i)] = rhob * sprim[velocityIdx(nmat,i)];
  }

  // consolidate primitives into state vector
  state.insert(state.end(), sprim.begin(), sprim.end());

  if (intsharp > 0 && srcFlag == 0)
  {
    std::vector< tk::real > vfmax(nmat, 0.0), vfmin(nmat, 0.0);

    tk::THINCReco(rdof, nmat, e, inpoel, coord, geoElem,
      ref_gp, U, P, intInd, matInt, vfmin, vfmax, state);
  }

  // physical constraints
  if (state.size() > ncomp) {
    for (std::size_t k=0; k<nmat; ++k) {
      state[ncomp+pressureIdx(nmat,k)] = constrain_pressure( mat_blk,
        state[ncomp+pressureIdx(nmat,k)], state[densityIdx(nmat,k)],
        state[volfracIdx(nmat,k)], k );
    }
  }

  return state;
}

void
safeReco( std::size_t rdof,
          std::size_t nmat,
          std::size_t el,
          int er,
          const Fields& U,
          std::array< std::vector< real >, 2 >& state )
// *****************************************************************************
//  Compute safe reconstructions near material interfaces
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nmat Total number of material is PDE system
//! \param[in] el Element on the left-side of face
//! \param[in] er Element on the right-side of face
//! \param[in] U Solution vector at recent time-stage
//! \param[in,out] state Second-order reconstructed state, at cell-face, that
//!   is being modified for safety
//! \details When the consistent limiting is applied, there is a possibility
//!   that the material densities and energies violate TVD bounds. This function
//!   enforces the TVD bounds locally
// *****************************************************************************
{
  using inciter::densityIdx;
  using inciter::energyIdx;
  using inciter::densityDofIdx;
  using inciter::energyDofIdx;

  if (er < 0) Throw("safe limiting cannot be called for boundary cells");

  auto eR = static_cast< std::size_t >(er);

  // define a lambda for the safe limiting
  auto safeLimit = [&]( std::size_t c, real ul, real ur )
  {
    // find min/max at the face
    auto uMin = std::min(ul, ur);
    auto uMax = std::max(ul, ur);

    // left-state limiting
    state[0][c] = std::min(uMax, std::max(uMin, state[0][c]));

    // right-state limiting
    state[1][c] = std::min(uMax, std::max(uMin, state[1][c]));
  };

  for (std::size_t k=0; k<nmat; ++k)
  {
    real ul(0.0), ur(0.0);

    // Density
    // establish left- and right-hand states
    ul = U(el, densityDofIdx(nmat, k, rdof, 0));
    ur = U(eR, densityDofIdx(nmat, k, rdof, 0));

    // limit reconstructed density
    safeLimit(densityIdx(nmat,k), ul, ur);

    // Energy
    // establish left- and right-hand states
    ul = U(el, energyDofIdx(nmat, k, rdof, 0));
    ur = U(eR, energyDofIdx(nmat, k, rdof, 0));

    // limit reconstructed energy
    safeLimit(energyIdx(nmat,k), ul, ur);
  }
}

} // tk::
