// *****************************************************************************
/*!
  \file      src/Inciter/AMR/Error.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Class for computing error estimates for mesh refinement
  \details   Class for computing error estimates for mesh refinement.
*/
// *****************************************************************************
#ifndef Error_h
#define Error_h

#include "Fields.hpp"
#include "Inciter/Options/AMRError.hpp"
#include "AMR/edge.hpp"

namespace AMR {

//! Class for computing error estimates for mesh refinement
class Error {

 private:
   using ncomp_t = tk::ncomp_t;

  public:
    //! Compute error estimate for a scalar quantity
    tk::real scalar( const tk::Fields& u,
                     const edge_t& edge,
                     ncomp_t c,
                     const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     const std::pair< std::vector< std::size_t >,
                                      std::vector< std::size_t > >& esup,
                     inciter::ctr::AMRErrorType err ) const;

  private:
    //! Estimate error for scalar quantity on edge based on jump in solution
    tk::real
    error_jump( const tk::Fields& u,
                const edge_t& edge,
                ncomp_t c ) const;

    //! Estimate error for scalar quantity on edge based on Hessian of solution
    tk::real
    error_hessian( const tk::Fields& u,
                   const edge_t& edge,
                   ncomp_t c,
                   const std::array< std::vector< tk::real >, 3 >& coord,
                   const std::vector< std::size_t >& inpoel,
                   const std::pair< std::vector< std::size_t >,
                                    std::vector< std::size_t > >& esup ) const;
};

} // AMR::

#endif // Error_h
