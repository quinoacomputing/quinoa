// *****************************************************************************
/*!
  \file      src/PDE/Reconstruction.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Reconstruction for reconstructed Galerkin methods
  \details   This file contains functions that reconstruct an "n"th order
    polynomial to an "n+1"th order polynomial using a least-squares
    reconstruction procedure, used for reconstructed discontinuous Galerkin (DG)
    methods. It also contains functions used to compute reconstruction in 1D,
    used in edge-based continuous Galerkin methods.
*/
// *****************************************************************************
#ifndef Reconstruction_h
#define Reconstruction_h

#include "Types.hpp"
#include "Fields.hpp"
#include "Limiter.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "FunctionPrototypes.hpp"
#include "Integrate/Basis.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "EoS/EOS.hpp"
#include "Kokkos_Core.hpp"
#include "KokkosDevice.hpp"

// Max size of quadrature array
// constexpr std::size_t NQUAD_MAX = 14;
// Max number of degrees of freedom
// constexpr std::size_t NDOF_MAX = 10;

// using execution_space = Kokkos::DefaultExecutionSpace;
// using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

namespace tk {

using ncomp_t = tk::ncomp_t;

//! \brief Reconstruct the second-order solution using least-squares approach
//!   from an extended stencil involving the node-neighbors
void
recoLeastSqExtStencil(
  std::size_t rdof,
  std::size_t e,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const Fields& geoElem,
  Fields& W,
  const std::vector< std::size_t >& varList );

//! Transform the reconstructed P1-derivatives to the Dubiner dofs
void
transform_P0P1( std::size_t rdof,
                std::size_t e,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                Fields& W,
                const std::vector< std::size_t >& varList );

//! Compute THINC reconstructions near material interfaces
void
THINCReco(
  std::size_t rdof,
  std::size_t nmat,
  std::size_t e,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const std::array< real, 3 >& xp,
  const Fields& U,
  const Fields& P,
  bool intInd,
  const std::vector< std::size_t >& matInt,
  const std::vector< real >& vfmin,
  const std::vector< real >& vfmax,
  std::vector< real >& state );

//! Compute THINC reconstructions for linear advection (transport)
void
THINCRecoTransport(
  std::size_t rdof,
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
  std::vector< real >& state );

//! Old THINC reconstruction function for volume fractions near interfaces
void
THINCFunction_old( std::size_t rdof,
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
  std::vector< real >& alReco );

KOKKOS_INLINE_FUNCTION 
void THINCFunction_old( std::size_t rdof,
                        std::size_t nmat,
                        std::size_t e,
                        std::size_t dof_e,
                        Kokkos::View<const size_t*, memory_space> inpoel,
                        Kokkos::View<const real*, memory_space> cx,
                        Kokkos::View<const real*, memory_space> cy,
                        Kokkos::View<const real*, memory_space> cz,
                        const Kokkos::Array<real, 3>& ref_xp,
                        real vol,
                        real bparam,
                        Kokkos::Array<real, NALSOL_MAX>& alSol,
                        bool intInd,
                        Kokkos::Array<size_t, NMAT_MAX>& matInt,
                        Kokkos::Array<real, NMAT_MAX>& alReco)
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
        bmod = bmod < (alk-epsl)/(epsh-epsl) * (bred - bparam) + bparam ? 
            bmod : (alk-epsl)/(epsh-epsl) * (bred - bparam) + bparam;
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

    // Extract the element coordinates
    Kokkos::Array<Kokkos::Array<real, 3>, 4> coordel;
    for (int i=0; i<4; i++) {
      coordel[i][0] = cx(inpoel(4*e + i));
      coordel[i][1] = cy(inpoel(4*e + i));
      coordel[i][2] = cz(inpoel(4*e + i));
    }

    Kokkos::Array<Kokkos::Array<real, NDOF_MAX>, 3> dBdx = {};
    Kokkos::Array<real, NDOF_MAX> B = {};
    auto jacInv =
      inverseJacobian(coordel[0], coordel[1], coordel[2], coordel[3] );

    eval_dBdx_p1(dof_e, jacInv, dBdx);
    
    Kokkos::Array<real, 3> nInt = {};
    Kokkos::Array<Kokkos::Array<real, 3>, NMAT_MAX> ref_n = {};

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
        Kokkos::Array<real, 3> axis = {
          coordel[i+1][0]-coordel[0][0],
          coordel[i+1][1]-coordel[0][1],
          coordel[i+1][2]-coordel[0][2]};
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
          Bc = std::exp(0.5*beta*(ref_n[k][1])+ref_n[k][2]);
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
        auto al_c = 0.5 * (1.0 + Kokkos::tanh(beta*(ref_n[k][0]*ref_xp[0]
          + ref_n[k][1]*ref_xp[1]+ref_n[k][2]*ref_xp[2] + d)));

        //! nested std::max, min might pose a problem
        auto m = min_lim > al_c ? min_lim : al_c;
        alReco[k] = max_lim < m ? max_lim : m;

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
};

KOKKOS_INLINE_FUNCTION void
THINCReco( std::size_t rdof,
           std::size_t nmat,
           std::size_t e, 
           std::size_t ncomp,
           std::size_t m_nprop,
           std::size_t p_nprop,
           std::size_t geo_nprop,
           std::size_t dof_e,
           tk::real bparam,
           Kokkos::View<const size_t*, memory_space> inpoel,
           Kokkos::View<const real*, memory_space> cx,
           Kokkos::View<const real*, memory_space> cy, 
           Kokkos::View<const real*, memory_space> cz,
           Kokkos::View<const real*, memory_space> geoElem,
           const Kokkos::Array<real, 3>& ref_xp,
           Kokkos::View<const real*, memory_space> U,
           Kokkos::View<const real*, memory_space> P,
           bool intInd,
           Kokkos::View<const size_t*, memory_space> solidx,
           Kokkos::Array<size_t, NMAT_MAX> matInt,
           [[maybe_unused]] const Kokkos::Array<real, NMAT_MAX>& vfmin,
           [[maybe_unused]] const Kokkos::Array<real, NMAT_MAX>& vfmax,
           Kokkos::Array<real, NSTATE_MAX>& state )
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

  // Step-1: Perform THINC reconstruction
  // create a vector of volume-fractions and pass it to the THINC function

  Kokkos::Array<real, NALSOL_MAX> alSol;
  Kokkos::Array<real, NMAT_MAX> alReco;
  for (std::size_t k=0; k<nmat; ++k) {
    auto mark = k*rdof;
    for (std::size_t i=0; i<rdof; ++i) {
      alSol[mark+i] = U(e * m_nprop + volfracDofIdx(nmat,k,rdof,i));
    }
    // initialize with TVD reconstructions which will be modified if near
    // material interface
    alReco[k] = state[volfracIdx(nmat,k)];
  }
  THINCFunction_old(rdof, nmat, e, dof_e, inpoel, cx, cy, cz, ref_xp, geoElem(e *geo_nprop), bparam,
    alSol, intInd, matInt, alReco);

  // check reconstructed volfracs for positivity
  bool neg_vf = false;
  for (std::size_t k=0; k<nmat; ++k) {
    if (alReco[k] < 1e-16 && matInt[k] > 0) neg_vf = true;
  }
  /*
  for (std::size_t k=0; k<nmat; ++k) {
    if (neg_vf) {
      printf("Material-id:        %d\n", k);
      printf("Volume-fraction:    %0.17f\n",alReco(k));
      printf("Cell-avg vol-frac:  %0.17f\n", U(e * m_nprop + volfracDofIdx(nmat,k,rdof,0)));
      printf("Material-interface? %d\n", intInd);
      printf("Mat-k-involved?     %d\n", matInt(k));
    }
  }
  */
  if (neg_vf) 
    printf("Material has negative volume fraction after THINC "
    "reconstruction.\n");

  // Step-2: Perform consistent reconstruction on other conserved quantities
  if (intInd)
  {
    auto rhobCC(0.0), rhobHO(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alCC = U(e * m_nprop + volfracDofIdx(nmat,k,rdof,0));
      alCC = 1e-14 > alCC ? 1e-14 : alCC;

      if (matInt[k])
      {
        state[volfracIdx(nmat,k)] = alReco[k];
        state[densityIdx(nmat,k)] = alReco[k]
          * U(e * m_nprop + densityDofIdx(nmat,k,rdof,0))/alCC;
        state[energyIdx(nmat,k)] = alReco[k]
          * U(e * m_nprop + energyDofIdx(nmat,k,rdof,0))/alCC;
        state[ncomp+pressureIdx(nmat,k)] = alReco[k]
          * P(e * p_nprop + pressureDofIdx(nmat,k,rdof,0))/alCC;
        if (solidx(k) > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              state[deformIdx(nmat,solidx(k),i,j)] =
                U(e * m_nprop + deformDofIdx(nmat,solidx(k),i,j,rdof,0));

          for (std::size_t i=0; i<6; ++i)
            state[ncomp+stressIdx(nmat,solidx(k),i)] = alReco[k]
              * P(e * p_nprop + stressDofIdx(nmat,solidx(k),i,rdof,0))/alCC;
        }
      }

      rhobCC += U(e * m_nprop + densityDofIdx(nmat,k,rdof,0));
      rhobHO += state[densityIdx(nmat,k)];
    }

    // consistent reconstruction for bulk momentum
    for (std::size_t i=0; i<3; ++i)
    {
      state[momentumIdx(nmat,i)] = rhobHO
        * U(e * m_nprop +  momentumDofIdx(nmat,i,rdof,0))/rhobCC;
      state[ncomp+velocityIdx(nmat,i)] =
        P(e * p_nprop +  velocityDofIdx(nmat,i,rdof,0));
    }
  }
}

//! New THINC reconstruction function for volume fractions near interfaces
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
  std::vector< real >& alReco );

//! Compute the temperatures based on FV conserved quantities
void
computeTemperaturesFV(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t nmat,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  const tk::Fields& unk,
  const tk::Fields& prim,
  const std::vector< int >& srcFlag,
  tk::Fields& T );

//! Evaluate polynomial solution at quadrature point
void
evalPolynomialSol(
  const std::vector< inciter::EOS >& mat_blk,
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
  const Fields& P,
  std::vector< tk::real >& state );

KOKKOS_INLINE_FUNCTION
// NOTE TO FUTURE INTERNS:
// I removed mat_blk from this device fn's param list
// mat_blk was never derefed here but capturing std::vector<EOS> into
// device lambda closure is not trivially copyable and costs host heap alloc
// per launch. If enforcePhysicalConstraints is reinstated/implemented then
// the EOS needs to be first mirrored into device-addressable POD
// That should be addressed by the union_EOS branch which has open pull req now
void evalPolynomialSol( //const std::vector< inciter::EOS >& mat_blk,
                        int intsharp,
                        std::size_t ncomp,
                        std::size_t nprim,
                        std::size_t rdof,
                        std::size_t nmat,
                        std::size_t e,
                        std::size_t dof_e,
                        std::size_t m_nprop,
                        std::size_t p_nprop,
                        std::size_t geo_nprop,
                        tk::real bparam,
                        Kokkos::View<const size_t*, memory_space> solidx,
                        Kokkos::View<const size_t*, memory_space> inpoel,
                        Kokkos::View<const tk::real*, memory_space> cx,
                        Kokkos::View<const tk::real*, memory_space> cy,
                        Kokkos::View<const tk::real*, memory_space> cz,
                        Kokkos::View<const tk::real*, memory_space> geoElem,
                        const Kokkos::Array<tk::real, 3>& ref_gp,
                        Kokkos::Array<tk::real, NDOF_MAX>& B,
                        Kokkos::View<const tk::real*, memory_space> U,
                        Kokkos::View<const tk::real*, memory_space> P,
                        Kokkos::Array<tk::real, NSTATE_MAX>& state )
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

  // since we are combining state, and primitives at the end,
  // just have one single "state" view and use subviews
  //auto state_sub = Kokkos::subview(state, Kokkos::make_pair(0, 1));
  //auto sprim_sub = Kokkos::subview(state, Kokkos::make_pair(1, 2));
  eval_state( ncomp, rdof, dof_e, e, m_nprop, U, B, state, 0); //?DONE
  eval_state( nprim, rdof, dof_e, e, p_nprop, P, B, state, ncomp); //?DONE

  // interface detection
  Kokkos::Array<size_t, NMAT_MAX> matInt = {};
  bool intInd(false);
  Kokkos::Array<real, NMAT_MAX> alAvg = {};

  if (nmat > 1) {
    for (std::size_t k=0; k<nmat; ++k) {
      alAvg[k] = U(e * m_nprop + inciter::volfracDofIdx(nmat,k,rdof,0)); //?DONE
      intInd = inciter::interfaceIndicator(nmat, alAvg, matInt); //?DONE
    }


  if (intsharp > 0)
  {
    Kokkos::Array<real, NMAT_MAX> vfmin = {};
    Kokkos::Array<real, NMAT_MAX> vfmax = {};

    // Until the appropriate setup for activating THINC with Transport
    // is ready, the following two chunks of code will need to be commented
    // for using THINC with Transport
    //for (std::size_t k=0; k<nmat; ++k) {
    //  vfmin[k] = VolFracMax(el, 2*k, 0);
    //  vfmax[k] = VolFracMax(el, 2*k+1, 0);
    //}
    tk::THINCReco(rdof, nmat, e, ncomp, m_nprop, p_nprop, geo_nprop, dof_e,
                  bparam, inpoel, cx, cy, cz, geoElem, ref_gp, U, P, intInd,
                  solidx, matInt, vfmin, vfmax, state);

    // Until the appropriate setup for activating THINC with Transport
    // is ready, the following lines will need to be uncommented for
    // using THINC with Transport
    //tk::THINCRecoTransport(rdof, nmat, el, inpoel, coord,
    //  geoElem, ref_gp_l, U, P, vfmin, vfmax, state[0]);
  }

  // // physical constraints
  // //? uncomment the below equations if pressure is an issue
  // //enforcePhysicalConstraints(mat_blk, ncomp, nmat, state);
  for (std::size_t k=0; k<nmat; ++k) {
      auto pres = state[ncomp+inciter::pressureIdx(nmat,k)];
      if (k == 0) {
        tk::real gamma = 4.4;
        tk::real pinf = 6.0e8;
        state[ncomp+inciter::pressureIdx(nmat,k)] = -gamma*pinf+1e-12 > pres ? -gamma*pinf+1e-12 : pres;
      }
      else
        state[ncomp+inciter::pressureIdx(nmat,k)] = 1e-12 > pres ? 1e-12 : pres;
    }
  }

}

//! Evaluate second-order FV solution at quadrature point
std::vector< tk::real >
evalFVSol(
  const std::vector< inciter::EOS >& mat_blk,
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
  int srcFlag );

//! Enforce physical constraints on state at quadrature point
void
enforcePhysicalConstraints(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t nmat,
  std::size_t ncomp,
  std::vector< tk::real >& state );

void
enforcePhysicalConstraints(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t nmat,
  std::size_t ncomp,
  Kokkos::View<real*, memory_space> state);

//! Evaluate mesh velocity at a quadrature point on a triangular face
std::array< real, 3 >
evaluateMeshVelTri(
  const std::size_t f,
  const std::size_t igp,
  const std::vector< std::size_t >& inpofa,
  const std::array< std::vector< real >, 2 >& coordgp,
  const Fields& W );

//! Evaluate mesh velocity at a quadrature point on a tetrahedral element
std::array< real, 3 >
evaluateMeshVelTet(
  const std::size_t e,
  const std::size_t igp,
  const std::vector< std::size_t >& inpoel,
  const std::array< std::vector< real >, 3 >& coordgp,
  const Fields& W );

//! Compute safe reconstructions near material interfaces
void
safeReco( std::size_t rdof,
          std::size_t nmat,
          std::size_t el,
          int er,
          const Fields& U,
          std::array< std::vector< real >, 2 >& state );

} // tk::

#endif // Reconstruction_h
