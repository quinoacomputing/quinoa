// *****************************************************************************
/*!
  \file      src/LoadBalance/LoadDistributor.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Load distributors
  \details   Load distributors compute chunksize based on the degree of
     virtualization.
*/
// *****************************************************************************

#include <limits>

#include "Types.h"
#include "LoadDistributor.h"
#include "Exception.h"

namespace tk {

uint64_t
linearLoadDistributor( tk::real virtualization,
                       uint64_t load,
                       int npe,
                       uint64_t& chunksize,
                       uint64_t& remainder )
// *****************************************************************************
//  Compute linear load distribution for given total work and virtualization
//! \param[in] virtualization Degree of virtualization [0.0...1.0]
//! \param[in] load Total load, e.g., number of particles, number of mesh cells
//! \param[in] npe Number of processing elements to distribute the load to
//! \param[inout] chunksize Chunk size, see detailed description
//! \param[inout] remainder Remainder, see detailed description
//! \return Number of work units
//! \details Compute load distibution (number of chares and chunksize) based on
//!   total work (e.g., total number of particles) and virtualization
//!
//!   The virtualization parameter, specified by the user, is a real number
//!   between 0.0 and 1.0, inclusive, which controls the degree of
//!   virtualization or over-decomposition. Independent of the value of
//!   virtualization the work is approximately evenly distributed among the
//!   available processing elements, given by npe. For zero virtualization (no
//!   over-decomposition), the work is simply decomposed into total_work/numPEs,
//!   which yields the smallest number of Charm++ chares and the largest chunks
//!   of work units. The other extreme is unity virtualization, which decomposes
//!   the total work into the smallest size work units possible, yielding the
//!   largest number of Charm++ chares. Obviously, the optimum will be between
//!   0.0 and 1.0, depending on the problem.
//!
//!   The formula implemented uses a linear relationship between the
//!   virtualization parameter and the number of work units with the extremes
//!   described above. The formula is given by
//!
//!   chunksize = (1 - n) * v + n;
//!
//!   where
//!    - v = degree of virtualization
//!    - n = load/npes
//!    - load = total work, e.g., number of particles, number of mesh cells
//!    - npes = number of hardware processing elements
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( virtualization > -std::numeric_limits< tk::real >::epsilon() &&
          virtualization < 1.0+std::numeric_limits< tk::real >::epsilon(),
          "Virtualization parameter must be between [0.0...1.0]" );
  Assert( npe > 0, "Number of processing elements must be larger than zero" );

  // Compute minimum number of work units
  const auto n = static_cast< tk::real >( load ) / npe;

  // Compute work unit size based on the linear formula above
  chunksize = static_cast< uint64_t >( (1.0 - n) * virtualization + n );

  // Compute number of work units with size computed ignoring remainder
  uint64_t nchare = load / chunksize;

  // Compute remainder of work if the above number of units were to be created
  remainder = load - nchare * chunksize;

  // Redistribute remainder among the work units for a more equal distribution
  chunksize += remainder / nchare;

  // Compute new remainder (after redistribution of the previous remainder)
  remainder = load - nchare * chunksize;

  // Return number of work units (number of Charm++ chares)
  return nchare;
}

} // tk::
