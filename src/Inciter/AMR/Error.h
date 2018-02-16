// *****************************************************************************
/*!
  \file      src/Inciter/AMR/Error.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Class for computing error estimates for mesh refinement
  \details   Class for computing error estimates for mesh refinement.
*/
// *****************************************************************************
#ifndef Error_h
#define Error_h

#include "Fields.h"
#include "Keywords.h"
#include "Inciter/Options/AMRError.h"

namespace AMR {

//! Class for computing error estimates for mesh refinement
class Error {

 private:
   using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Compute error estimate for a scalar quantity
    tk::real scalar( const tk::Fields& u,
                     const std::pair< std::size_t, std::size_t >& edge,
                     ncomp_t c,
                     const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     const std::pair< std::vector< std::size_t >,
                                      std::vector< std::size_t > >& esup,
                     inciter::ctr::AMRErrorType err );

  private:
    //! Estimate error for scalar quantity on edge based on jump in solution
    tk::real
    error_jump( const tk::Fields& u,
                const std::pair< std::size_t, std::size_t >& edge,
                ncomp_t c );

    //! Estimate error for scalar quantity on edge based on Hessian of solution
    tk::real
    error_hessian( const tk::Fields& u,
                   const std::pair< std::size_t, std::size_t >& edge,
                   ncomp_t c,
                   const std::array< std::vector< tk::real >, 3 >& coord,
                   const std::vector< std::size_t >& inpoel,
                   const std::pair< std::vector< std::size_t >,
                                    std::vector< std::size_t > >& esup );
};

} // AMR::

#endif // Error_h
