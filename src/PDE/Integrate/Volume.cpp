// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing volume integrals for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing volume integrals for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************

#include "Volume.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "Reconstruction.hpp"
#include "MultiMatTerms.hpp"
#include "Kokkos_Core.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

using execution_space = Kokkos::Serial;
using memory_space = Kokkos::HostSpace;
using range_policy = Kokkos::RangePolicy<execution_space>;
using UnManagedMem =Kokkos::MemoryTraits<Kokkos::Unmanaged>;

template <typename T>
auto changeToView(T* object, size_t n) {
    Kokkos::View<T*, Kokkos::HostSpace, UnManagedMem> object_view(object, n);
    return object_view;
}

namespace inciter {
  extern ctr::InputDeck g_inputdeck; //! is this allowed?
};

void
tk::volInt( std::size_t nmat,
            real t,
            const std::vector< inciter::EOS >& mat_blk,
            const std::size_t ndof,
            const std::size_t rdof,
            const std::size_t nelem,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            const Fields& geoElem,
            const FluxFn& flux,
            const VelFn& vel,
            const Fields& U,
            const Fields& P,
            const std::vector< std::size_t >& ndofel,
            Fields& R,
            int intsharp )
// *****************************************************************************
//  Compute volume integrals for DG
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] t Physical time
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Total number of degrees of freedom included reconstructed ones
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{

  int ncomp = U.nprop()/rdof;
  int nprim = P.nprop()/rdof;

  size_t m_nprop = U.nprop();

  // solidx is a vector
  const auto& solidx = inciter::g_inputdeck.get<
      tag::matidxmap, tag::solidx >();

  // bparm is a scalar
  auto bparam = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp_param >();

  // compute volume integrals
  Kokkos::initialize();
  {
    //! Transfer all the constants to View is not necessary!
    //Transfer solidx is a vector
    auto solidx_h_view = changeToView(solidx.data(), nmat);
    Kokkos::View<size_t*, memory_space> solidx_d_view("solidx", nmat);
    Kokkos::deep_copy(solidx_d_view, solidx_h_view);

    // Transfer inpoel variable
    size_t inpoel_size = inpoel.size();
    auto inpoel_h_view = changeToView(inpoel.data(), inpoel_size);
    Kokkos::View<size_t*, memory_space> inpoel_d_view("inpoel device view", inpoel_size);
    Kokkos::deep_copy(inpoel_d_view, inpoel_h_view);

    //Transfer coord (nodal coordinates)

    size_t coordx_size = coord[0].size();
    auto cx_h_view = changeToView(coord[0].data(), coordx_size);
    Kokkos::View<real*, memory_space> cx_d_view("cx device view", coordx_size);
    Kokkos::deep_copy(cx_d_view, cx_h_view);

    size_t coordy_size = coord[1].size();
    auto cy_h_view = changeToView(coord[1].data(), coordy_size);
    Kokkos::View<real*, memory_space> cy_d_view("cy device view", coordy_size);
    Kokkos::deep_copy(cy_d_view, cy_h_view);

    size_t coordz_size = coord[2].size();
    auto cz_h_view = changeToView(coord[2].data(), coordz_size);
    Kokkos::View<real*, memory_space> cz_d_view("cz device view", coordz_size);
    Kokkos::deep_copy(cz_d_view, cz_h_view);

    // nodefel transfer

    size_t ndofel_size = ndofel.size();
    auto ndofel_h_view = changeToView(ndofel.data(), ndofel_size);
    Kokkos::View<double*, memory_space> ndofel_d_view("nodefel device view", ndofel_size);
    Kokkos::deep_copy(ndofel_d_view, ndofel_h_view);

    // geoElem, U, P, R transfer

    size_t geoElem_size = geoElem.getSize();
    Kokkos::View<real*, memory_space> geoElem_d_view("geoElem_d_view", geoElem_size);
    auto geoElem_h_view = changeToView(geoElem.getPointer(), geoElem_size);
    Kokkos::deep_copy(geoElem_d_view, geoElem_h_view);

    size_t P_size = P.getSize();
    Kokkos::View<real*, memory_space> P_d_view("P_d_view", P_size);
    auto P_h_view = changeToView(P.getPointer(), P_size);
    Kokkos::deep_copy(P_d_view, P_h_view);

    size_t U_size = U.getSize();
    Kokkos::View<real*, memory_space> U_d_view("U_d_view", U_size);
    auto U_h_view = changeToView(U.getPointer(), U_size);
    Kokkos::deep_copy(U_d_view, U_h_view);

    size_t R_size = R.getSize();
    Kokkos::View<real*, memory_space> R_d_view("R_d_view", R_size);
    auto R_h_view = changeToView(R.getPointer(), R_size);
    Kokkos::deep_copy(R_d_view, R_h_view);

    //create View variables in device space that will be used inside the kernel (parallel env)
   // The sizes are only known inside the kernel!!
    Kokkos::View<real**, memory_space> coordgp("coordgp_d_view", 3, 3);
    Kokkos::View<real*, memory_space> wgp("wgp_d_view", 1);
    Kokkos::View<real**, memory_space> dBdx("dBdx_d_view", 3, 3);
    Kokkos::View<real*, memory_space> B("B", 3);

    // for flux evaluation, but size is known before kernel
    Kokkos::View<real***, memory_space> g("g_d_view", nmat, 3, 3);
    Kokkos::View<real***, memory_space> asig("asig_d_view", nmat, 3, 3);
    Kokkos::View<real*, memory_space> al("al", nmat);
    Kokkos::View<real**, memory_space> fl("fl", ncomp, 3);
    Kokkos::View<real*, memory_space> apk("apk", nmat);

    //Need for evalPolynomialSol function
    Kokkos::View<real*, memory_space> state("state", 2*ncomp); // state has state + sprim length
    Kokkos::View<size_t*, memory_space> matInt("matInt", nmat);
    Kokkos::View<real*, memory_space> alAvg("alAvg", nmat);
    Kokkos::View<real*, memory_space> vfmax("vfmax", nmat);
    Kokkos::View<real*, memory_space> vfmin("vfmin", nmat);

    //Need for THINC, but also evalPolynomialSol
    Kokkos::View<real*, memory_space> alSol("alSol", rdof*nmat);
    Kokkos::View<real*, memory_space> alReco("alReco", nmat);
    Kokkos::View<Kokkos::Array<real, 3>*, memory_space> ref_n("ref_n", nmat);

    Kokkos::parallel_for(range_policy(0, nelem), KOKKOS_LAMBDA(const size_t e)
    {
        if(ndofel_d_view(e) > 1)
        {
          auto ng = tk::NGvol(ndofel_d_view(e)); //?DONE

          // arrays for quadrature points /
          //! cant' resize in paralle for
          //!Kokkos::resize(coordgp, 3, static_cast<const std::size_t>(ng));
          //!Kokkos::resize(wgp, static_cast<const std::size_t>(ng));

          GaussQuadratureTet(ng, coordgp, wgp ); //?DONE

          // Extract the element coordinates
          
          /*std::array< std::array< real, 3>, 4 > coordel {{
            {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
            {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
            {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
            {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
          }};
          */

          Kokkos::Array<Kokkos::Array<real, 3>, 4> coordel;
          for (int i = 0;i < 4;i++) {
              coordel[i][0] = cx_d_view(inpoel_d_view(4*e + i));
              coordel[i][1] = cy_d_view(inpoel_d_view(4*e + i));
              coordel[i][2] = cz_d_view(inpoel_d_view(4*e + i));
            }
        
          //jacInv is Kokkos::Array based matrix
          auto jacInv =
                  inverseJacobian(coordel[0], coordel[1], coordel[2], coordel[3] ); // ?DONE

          auto dof_el = ndofel_d_view(e); //? DONE

          //! Pass in dBdx rather than returning it since I have already created dBdx as view type
          eval_dBdx_p1(dof_el, jacInv, dBdx); //?DONE
          
            // Gaussian quadrature
            for (std::size_t igp=0; igp<ng; ++igp)
            {
              if (dof_el > 4)
                eval_dBdx_p2( igp, coordgp, jacInv, dBdx); //?DONE

              // Compute the coordinates of quadrature point at physical domain
              auto gp = eval_gp( igp, coordel, coordgp); // ?DONE

              // Compute the basis function
              // B is Kokkos::View //TODO: Fix B
              eval_basis( dof_el, coordgp(0, igp), coordgp(1, igp),
                                  coordgp(2, igp), B); //?DONE

              auto wt = wgp(igp) * geoElem_d_view(e * m_nprop);  //?DONE

              evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
                rdof, nmat, e, ndofel_d_view(e), m_nprop, bparam, solidx_d_view, inpoel_d_view, cx_d_view, cy_d_view, cz_d_view, geoElem_d_view, // * pass in cx, cy, cx rather
                {{coordgp(0, igp), coordgp(1, igp), coordgp(2, igp)}}, B, U_d_view, P_d_view, state,
                matInt, alAvg, vfmax, vfmin, alSol, alReco, dBdx, ref_n);

              // evaluate prescribed velocity (if any)
              //auto v = vel( ncomp, gp[0], gp[1], gp[2], t ); 

              // comput flux
              
              tk::fluxTerms_multimat_kokkos(ncomp, nmat, solidx_d_view, 
                  mat_blk, state, g, asig, al, fl, apk);

              update_rhs(ncomp, ndof, dof_el, wt, m_nprop, e, dBdx, fl, R_d_view); //?DONE
          }
        }
      });
    };
     Kokkos::finalize();
}
 

void tk::update_rhs( ncomp_t ncomp,
                const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt, 
                const std::size_t e,
                const std::array< std::vector<tk::real>, 3 >& dBdx,
                const std::vector< std::array< tk::real, 3 > >& fl,
                Fields& R )
// *****************************************************************************
//  Update the rhs by adding the source term integrals
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] dBdx Vector of basis function derivatives
//! \param[in] fl Vector of numerical flux
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( dBdx[0].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( dBdx[1].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( dBdx[2].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( fl.size() == ncomp, "Size mismatch for flux term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(e, mark+1) +=
      wt * (fl[c][0]*dBdx[0][1] + fl[c][1]*dBdx[1][1] + fl[c][2]*dBdx[2][1]);
    R(e, mark+2) +=
      wt * (fl[c][0]*dBdx[0][2] + fl[c][1]*dBdx[1][2] + fl[c][2]*dBdx[2][2]);
    R(e, mark+3) +=
      wt * (fl[c][0]*dBdx[0][3] + fl[c][1]*dBdx[1][3] + fl[c][2]*dBdx[2][3]);

    if( ndof_el > 4 )
    {
      R(e, mark+4) +=
        wt * (fl[c][0]*dBdx[0][4] + fl[c][1]*dBdx[1][4] + fl[c][2]*dBdx[2][4]);
      R(e, mark+5) +=
        wt * (fl[c][0]*dBdx[0][5] + fl[c][1]*dBdx[1][5] + fl[c][2]*dBdx[2][5]);
      R(e, mark+6) +=
        wt * (fl[c][0]*dBdx[0][6] + fl[c][1]*dBdx[1][6] + fl[c][2]*dBdx[2][6]);
      R(e, mark+7) +=
        wt * (fl[c][0]*dBdx[0][7] + fl[c][1]*dBdx[1][7] + fl[c][2]*dBdx[2][7]);
      R(e, mark+8) +=
        wt * (fl[c][0]*dBdx[0][8] + fl[c][1]*dBdx[1][8] + fl[c][2]*dBdx[2][8]);
      R(e, mark+9) +=
        wt * (fl[c][0]*dBdx[0][9] + fl[c][1]*dBdx[1][9] + fl[c][2]*dBdx[2][9]);
    }
  }
}


//! overloaded version of update_rhs for Kokkos
KOKKOS_INLINE_FUNCTION
void tk::update_rhs( ncomp_t ncomp,
                const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt,
                const std::size_t m_nprop,
                const std::size_t e,
                Kokkos::View<const real**, memory_space> dBdx,
                Kokkos::View<const real**, memory_space> fl,
                Kokkos::View<real*, memory_space> R)
 {
  for (ncomp_t c=0; c<ncomp; ++c)
    {
      auto mark = c*ndof;
      R(e * m_nprop + mark+1) +=
        wt * (fl(c, 0)*dBdx(0, 1) + fl(c, 1)*dBdx(1, 1) + fl(c, 2)*dBdx(2, 1));
      R(e * m_nprop + mark+2) +=
        wt * (fl(c, 0)*dBdx(0, 2) + fl(c, 1)*dBdx(1, 2) + fl(c, 2)*dBdx(2, 2));
      R(e * m_nprop + mark+3) +=
        wt * (fl(c, 0)*dBdx(0, 3) + fl(c, 1)*dBdx(1, 3) + fl(c, 2)*dBdx(2, 3));

      if( ndof_el > 4 )
      {
        R(e * m_nprop + mark+4) +=
          wt * (fl(c, 0)*dBdx(0, 4) + fl(c, 1)*dBdx(1, 4) + fl(c, 2)*dBdx(2, 4));
        R(e * m_nprop + mark+5) +=
          wt * (fl(c, 0)*dBdx(0, 5) + fl(c, 1)*dBdx(1, 5) + fl(c, 2)*dBdx(2, 5));
        R(e * m_nprop + mark+6) +=
          wt * (fl(c, 0)*dBdx(0, 6) + fl(c, 1)*dBdx(1, 6) + fl(c, 2)*dBdx(2, 6));
        R(e * m_nprop + mark+7) +=
          wt * (fl(c, 0)*dBdx(0, 7) + fl(c, 1)*dBdx(1, 7) + fl(c, 2)*dBdx(2, 7));
        R(e * m_nprop + mark+8) +=
          wt * (fl(c, 0)*dBdx(0, 8) + fl(c, 1)*dBdx(1, 8) + fl(c, 2)*dBdx(2, 8));
        R(e * m_nprop + mark+9) +=
          wt * (fl(c, 0)*dBdx(0, 9) + fl(c, 1)*dBdx(1, 9) + fl(c, 2)*dBdx(2, 9));
      }
    }
  }
