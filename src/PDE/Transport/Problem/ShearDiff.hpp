// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/ShearDiff.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for scalar transport equations
  \details   This file declares a Problem policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) and PDE/Transport/DGTransport.h
    implementing cell-centered discontinuous Galerkin (DG) discretizations.
    See PDE/Transport/Problem.h for general requirements on Problem policy
    classes for cg::Transport and dg::Transport.
*/
// *****************************************************************************
#ifndef TransportProblemShearDiff_h
#define TransportProblemShearDiff_h

#include <vector>
#include <array>

#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "Inciter/Options/Problem.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

/*! Transport PDE problem: diffusion of a shear layer
    \details This class implements the analytical solutions for the test
    problem, adopted from Okubo Akira Karweit Michael J. , (1969),
    [DIFFUSION FROM A CONTINUOUS SOURCE IN A UNIFORM SHEAR
    FLOW](http://onlinelibrary.wiley.com/doi/10.4319/lo.1969.14.4.0514/abstract),
    Limnology and Oceanography, 14, doi: 10.4319/lo.1969.14.4.0514. In essence,
    this is a test problem for the advection-diffusion equation in 3D where the
    analytical solution is known in a closed form as the solution evolves in
    time. The initial solution is a Gaussian that is advected and diffused in
    time with an imposed constant-in-time velocity field that features
    advection and shear. Also, the diffusion coefficients can be different in
    the three coordinate directions. Note that t0 as well as all three
    components of the diffusion must be larger than zero at t=t0 to have a
    well-defined initial condition.
   
    In a nutshell, the equation solved is
    \f[
      \frac{\partial S}{\partial t} + \left(V_0 + \Omega_y y + \Omega_z z
      \right) \frac{\partial S}{\partial x} =
      A_x \frac{\partial^2S}{\partial x^2} +
      A_y \frac{\partial^2S}{\partial y^2} +
      A_z \frac{\partial^2S}{\partial z^2}
    \f]
    whose solution is given by
    \f[
      S(t,x,y,z,) = \frac{1}{8\pi^{3/2}(A_xA_yA_z)^{1/2}t^{3/2}
                             (1+\phi_3^2t^2)^{1/2}}
                    \exp\left[ -\frac{x-V_0t-(\Omega_yy+\Omega_zz)^2/2}
                                     {4A_xt(1+\phi_3^2t^2}
                               -\frac{y^2}{4A_yt}
                               -\frac{z^2}{4A_zt} \right]
    \f]
    where \f$ \phi_3^2 = (\Omega_y^2A_y/A_x + \Omega_z^2A_z/A_x)/12\f$.
    See also the paper.
*/
class TransportProblemShearDiff {
  private:
    using ncomp_t = tk::ncomp_t;
    using eq = newtag::transport;

  public:
    //! Initialize numerical solution
    static std::vector< tk::real >
    initialize( ncomp_t ncomp,
                const std::vector< EOS >& mat_blk, tk::real x, tk::real y,
                tk::real z, tk::real t );

    //! Evaluate analytical solution at (x,y,z,t) for all components
    static std::vector< tk::real >
    analyticSolution( ncomp_t ncomp,
                      const std::vector< EOS >& mat_blk, tk::real x,
                      tk::real y, tk::real z, tk::real t )
    { return initialize( ncomp, mat_blk, x, y, z, t ); }

    //! Do error checking on PDE parameters
    void errchk( ncomp_t ncomp ) const;

    //! Assign prescribed shear velocity at a point
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( ncomp_t ncomp,
                        tk::real,
                        tk::real y,
                        tk::real z,
                        tk::real );

    //! Return problem type
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SHEAR_DIFF; }
};

} // inciter::

#endif // TransportProblemShearDiff_h
