// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/CGAdvDiff.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for advection-diffusion using continuous Galerkin
  \details   This file declares a Physics policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) discretizations.
    See PDE/Transport/Physics/CG.h for general requirements on Physics policy
    classes for cg::Transport.
*/
// *****************************************************************************
#ifndef TransportPhysicsCGAdvDiff_h
#define TransportPhysicsCGAdvDiff_h

#include "Types.hpp"
#include "Fields.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Physics.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace cg {

//! Physics policy for advection-diffusion using continuous Galerkin
class TransportPhysicsAdvDiff {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::transport;

  public:
    //! Add diffusion contribution to rhs
    //! \tparam Op Operation to specify boundary vs internal nodes contribution    
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this PDE
    //! \param[in] J Element Jacobi determinant
    //! \param[in] gid Local->global node id map
    //! \param[in] bid Local chare-boundary node ids (value) associated to global
    //!   node ids (key)
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] N Element node indices
    //! \param[in] u Solution at element nodes at recent time step
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    template< class Op >
    void
    diffusionRhs( ncomp_t e,
                  ncomp_t ncomp,
                  tk::real J,
                  const std::vector< std::size_t >& gid,
                  const std::unordered_map< std::size_t, std::size_t >& bid,
                  const std::array< std::array< tk::real, 3 >, 4 >& grad,
                  const std::array< std::size_t, 4 >& N,
                  const std::vector< std::array< tk::real, 4 > >& u,
                  const std::vector< const tk::real* >& r,
                  tk::Fields& R,
                  Op op ) const
    {
      // diffusivities for all components
      const auto& diff =
        g_inputdeck.get< tag::param, eq, tag::diffusivity >()[e];
    
      // add diffusion contribution to right hand side
      const auto d = J/6.0;
      for (std::size_t a=0; a<4; ++a) {
        if (op( bid.find(gid[N[a]]), end(bid) ))
          for (ncomp_t c=0; c<ncomp; ++c)
            for (std::size_t k=0; k<3; ++k) {
              const auto D = diff[ 3*c+k ];
                for (std::size_t b=0; b<4; ++b)
                  R.var(r[c],N[a]) -= d * D * grad[a][k] * grad[b][k] * u[c][b];
            }
      }
    }

    //! Compute the minimum time step size based on the diffusion
    tk::real
    diffusion_dt( ncomp_t e,
                  ncomp_t ncomp,
                  tk::real L,
                  const std::vector< std::array< tk::real, 4 > >& ) const;

    //! Return physics type
    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVDIFF; }
};

} // cg::
} // inciter::

#endif // TransportPhysicsCGAdvDiff_h
