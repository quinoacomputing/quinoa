// *****************************************************************************
/*!
  \file      src/PDE/Riemann/Upwind.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Upwind Riemann flux function
  \details   This file implements the upwind Riemann solver.
*/
// *****************************************************************************
#ifndef Upwind_h
#define Upwind_h

#include <vector>

#include "Types.hpp"
#include "Fields.hpp"
#include "Vector.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EosVariant.hpp"

namespace inciter {

//! Upwind Riemann solver
struct Upwind {

    //! Upwind Riemann solver flux function
    //! \param[in] fn Face/Surface normal
    //! \param[in] u Left and right unknown/state vector
    //! \param[in] v Prescribed velocity evaluated at the integration point
    //!   where this flux function is used for all scalar components in the
    //!   system of PDEs integrated
    //! \return Riemann solution using a central difference method
    //! \note The function signature must follow tk::RiemannFluxFn
    static tk::RiemannFluxFn::result_type
    flux( const std::vector< EOS >,
          const std::array< tk::real, 3 >& fn,
          const std::array< std::vector< tk::real >, 2 >& u,
          const std::vector< std::array< tk::real, 3 > >& v )
    {
      if (u[0].size() != v.size()) Throw( "Unequal components in transport" );
      std::vector< tk::real > flx( u[0].size(), 0 );

      for(std::size_t c=0; c<v.size(); ++c)
      {
        // wave speed based on prescribed velocity
        auto swave = tk::dot( v[c], fn );
      
        // upwinding
        tk::real splus  = 0.5 * (swave + fabs(swave));
        tk::real sminus = 0.5 * (swave - fabs(swave));
      
        flx[c] = splus * u[0][c] + sminus * u[1][c];
      }
      
      return flx;
    }
  
    //! Flux type accessor
    //! \return Flux type
    static ctr::FluxType type() noexcept { return ctr::FluxType::UPWIND; }
};

} // inciter::

#endif // Upwind_h
