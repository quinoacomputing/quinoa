// *****************************************************************************
/*!
  \file      src/Particles/Tracker.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Tracker tracks Lagrangian particles in physical space
  \details   Tracker tracks Lagrangian particles in physical space. It works on
    a chunk of the Eulerian mesh, and tracks particles in elements and across
    mesh chunks held by different Charm++ chares.
*/
// *****************************************************************************

#include "NoWarning/threefry.h"

#include "Random123.h"

#include "Tracker.h"
#include "Exception.h"

using tk::Tracker;

void
Tracker::genpar( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 std::size_t nchare,
                 int chid )
// *****************************************************************************
//  Generate particles to each of our mesh cells
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] nchare Total number of holder array chares
//! \param[in] chid Host chare ID (thisIndex)
//! \author F.J. Gonzalez
// *****************************************************************************
{
  Assert( m_elp.size() >= m_particles.nunk(),
          "Element-of-particle array not large enough" );

  auto rng = tk::Random123< r123::Threefry2x64 >( nchare );

  // Create a reference of mesh point coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Generate npar number of particles into each mesh cell
  auto npar = m_particles.nunk() / (inpoel.size()/4);
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    for (std::size_t p=0; p<npar; ++p) {
      std::array< tk::real, 4 > N;
      rng.uniform( chid, 3, N.data() );
      N[3] = 1.0 - N[0] - N[1] - N[2];
      if ( std::min(N[0],1.0-N[0]) > 0.0 && std::min(N[1],1.0-N[1]) > 0.0 &&
           std::min(N[2],1.0-N[2]) > 0.0 && std::min(N[3],1.0-N[3]) > 0.0 ) {
        const auto A = inpoel[e*4+0];
        const auto B = inpoel[e*4+1];
        const auto C = inpoel[e*4+2];
        const auto D = inpoel[e*4+3];
        const auto i = e * npar + p;
        m_particles(i,0,0) = x[A]*N[0] + x[B]*N[1] + x[C]*N[2] + x[D]*N[3];
        m_particles(i,1,0) = y[A]*N[0] + y[B]*N[1] + y[C]*N[2] + y[D]*N[3];
        m_particles(i,2,0) = z[A]*N[0] + z[B]*N[1] + z[C]*N[2] + z[D]*N[3];
        m_elp[i] = e;
      } else --p; // retry if particle was not generated into cell
    }
  }
}

std::vector< std::size_t >
Tracker::addpar( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const std::vector< std::size_t >& miss,
                 const std::vector< std::vector< tk::real > >& ps )
// *****************************************************************************
//  Try to find particles and add those found to the list of ours
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] miss Indices of particles to find
//! \param[in] ps Particle data associated to those particle indices to find
//! \return Particle indices found
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( ps.size() == miss.size(), "Size mismatch" );

  std::vector< std::size_t > found; // will store indices of particles found

  // try to find particles received
  for (std::size_t i=0; i<ps.size(); ++i)
    for (std::size_t e=0; e<inpoel.size()/4; ++e) {
      std::array< tk::real, 4 > N;
      auto last = m_particles.nunk();
      m_particles.push_back( ps[i] );
      if (parinel( coord, inpoel, last, e, N )) {
        found.push_back( miss[i] );
        m_elp.resize( last+1 );
        m_elp[ last ] = e;
        e = inpoel.size()/4;
      } else {
        m_particles.rm( { last } );
      }
    }

  return found;
}

bool
Tracker::parinel( const std::array< std::vector< tk::real >, 3 >& coord,
                  const std::vector< std::size_t >& inpoel,
                  std::size_t p,
                  std::size_t e,
                  std::array< tk::real, 4 >& N )
// *****************************************************************************
//  Search particle in a single mesh cell
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] p Particle index
//! \param[in] e Mesh cell index
//! \param[in,out] N Shapefunctions evaluated at the particle position
//! \return True if particle is in mesh cell
//! \author F.J. Gonzalez
// *****************************************************************************
{
  // Tetrahedron node indices
  const auto A = inpoel[e*4+0];
  const auto B = inpoel[e*4+1];
  const auto C = inpoel[e*4+2]; 
  const auto D = inpoel[e*4+3];

  // Tetrahedron node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Particle coordinates
  const auto& xp = m_particles(p,0,0);
  const auto& yp = m_particles(p,1,0);
  const auto& zp = m_particles(p,2,0);

  // Evaluate linear shapefunctions at particle locations using Cramer's Rule
  //    | xp |   | x1 x2 x3 x4 |   | N1 |
  //    | yp | = | y1 y2 y3 y4 | • | N2 |
  //    | zp |   | z1 z2 z3 z4 |   | N3 |
  //    | 1  |   | 1  1  1  1  |   | N4 |

  tk::real DetX = (y[B]*z[C] - y[C]*z[B] - y[B]*z[D] + y[D]*z[B] + 
    y[C]*z[D] - y[D]*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
    x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
    x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[C]*y[A]*z[D] + x[C]*y[D]*z[A] +
    x[D]*y[A]*z[C] - x[D]*y[C]*z[A] - x[B]*y[C]*z[D] + x[B]*y[D]*z[C] +
    x[C]*y[B]*z[D] - x[C]*y[D]*z[B] - x[D]*y[B]*z[C] + x[D]*y[C]*z[B];

  tk::real DetX1 = (y[D]*z[C] - y[C]*z[D] + y[C]*zp - yp*z[C] -
    y[D]*zp + yp*z[D])*x[B] + x[C]*y[B]*z[D] - x[C]*y[D]*z[B] - 
    x[D]*y[B]*z[C] + x[D]*y[C]*z[B] - x[C]*y[B]*zp + x[C]*yp*z[B] +
    xp*y[B]*z[C] - xp*y[C]*z[B] + x[D]*y[B]*zp - x[D]*yp*z[B] - 
    xp*y[B]*z[D] + xp*y[D]*z[B] + x[C]*y[D]*zp - x[C]*yp*z[D] - 
    x[D]*y[C]*zp + x[D]*yp*z[C] + xp*y[C]*z[D] - xp*y[D]*z[C];

  tk::real DetX2 = (y[C]*z[D] - y[D]*z[C] - y[C]*zp + yp*z[C] +
    y[D]*zp - yp*z[D])*x[A] + x[C]*y[D]*z[A] - x[C]*y[A]*z[D] +
    x[D]*y[A]*z[C] - x[D]*y[C]*z[A] + x[C]*y[A]*zp - x[C]*yp*z[A] -
    xp*y[A]*z[C] + xp*y[C]*z[A] - x[D]*y[A]*zp + x[D]*yp*z[A] +
    xp*y[A]*z[D] - xp*y[D]*z[A] - x[C]*y[D]*zp + x[C]*yp*z[D] + 
    x[D]*y[C]*zp - x[D]*yp*z[C] - xp*y[C]*z[D] + xp*y[D]*z[C];

  tk::real DetX3 = (y[D]*z[B] - y[B]*z[D] + y[B]*zp - yp*z[B] -
    y[D]*zp + yp*z[D])*x[A] + x[B]*y[A]*z[D] - x[B]*y[D]*z[A] -
    x[D]*y[A]*z[B] + x[D]*y[B]*z[A] - x[B]*y[A]*zp + x[B]*yp*z[A] +
    xp*y[A]*z[B] - xp*y[B]*z[A] + x[D]*y[A]*zp - x[D]*yp*z[A] - 
    xp*y[A]*z[D] + xp*y[D]*z[A] + x[B]*y[D]*zp - x[B]*yp*z[D] - 
    x[D]*y[B]*zp + x[D]*yp*z[B] + xp*y[B]*z[D] - xp*y[D]*z[B];

  tk::real DetX4 = (y[B]*z[C] - y[C]*z[B] - y[B]*zp + yp*z[B] +
    y[C]*zp - yp*z[C])*x[A] + x[B]*y[C]*z[A] - x[B]*y[A]*z[C] +
    x[C]*y[A]*z[B] - x[C]*y[B]*z[A] + x[B]*y[A]*zp - x[B]*yp*z[A] -
    xp*y[A]*z[B] + xp*y[B]*z[A] - x[C]*y[A]*zp + x[C]*yp*z[A] +
    xp*y[A]*z[C] - xp*y[C]*z[A] - x[B]*y[C]*zp + x[B]*yp*z[C] +
    x[C]*y[B]*zp - x[C]*yp*z[B] - xp*y[B]*z[C] + xp*y[C]*z[B];

  // Shape functions evaluated at particle location
  N[0] = DetX1/DetX;
  N[1] = DetX2/DetX;
  N[2] = DetX3/DetX;
  N[3] = DetX4/DetX;

  // if min( N^i, 1-N^i ) > 0 for all i, particle is in cell
  if ( std::min(N[0],1.0-N[0]) > 0 && std::min(N[1],1.0-N[1]) > 0 &&
       std::min(N[2],1.0-N[2]) > 0 && std::min(N[3],1.0-N[3]) > 0 )
  {
    m_elp.resize( p+1 );
    m_elp[ p ] = e; // store element of particle
    return true;
  } else {
    return false;
  }
}

void
Tracker::applyParBC( std::size_t i )
// *****************************************************************************
// Apply boundary conditions to particles
//! \author F.J. Gonzalez
// *****************************************************************************
{
  auto& x = m_particles(i,0,0);
  auto& y = m_particles(i,1,0);
  auto& z = m_particles(i,2,0);

  if (z > 1.0) z = 0.99;
  if (z < 0.0) z = 0.01;
  if (y > 1.0) y = 0.99;
  if (y < 0.0) y = 0.01;
  if (x > 1.0) x = 0.99;
  if (x < 0.0) x = 0.01;
}

void
Tracker::remove( const std::set< std::size_t >& idx )
// *****************************************************************************
// Remove particles
//! \param[in] idx Set of particle indices whose data to remove
//! \author J. Bakosi
// *****************************************************************************
{
  m_particles.rm( idx );

  auto rem = [ &idx ]( std::size_t i ) -> bool {
    if (idx.find(i) != end(idx)) return true;
    return false;
  };

  std::size_t last = 0;
  for(std::size_t i=0; i<m_elp.size(); ++i, ++last) {
    while( rem(i) ) ++i;
    if (i >= m_elp.size()) break;
    m_elp[last] = m_elp[i];
  }
  m_elp.resize( last );

  Assert( m_particles.nunk() == m_elp.size(),
          "Number of particles and the number of host elements unequal" );
}
