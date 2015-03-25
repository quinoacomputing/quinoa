//******************************************************************************
/*!
  \file      src/Mesh/DerivedData.C
  \author    J. Bakosi
  \date      Wed 25 Mar 2015 12:46:44 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
//******************************************************************************

#include <algorithm>

#include <make_unique.h>
#include <DerivedData.h>

namespace tk {

void
shiftToZero( std::vector< int >& inpoel )
//******************************************************************************
//  Shift node IDs to start with zero in element connectivity
//! \param[inout] inpoel Inteconnectivity of points and elements
//! \note It is okay to call this function with an empty container; it will
//!    simply return without throwing an exception.
//! \author J. Bakosi
//******************************************************************************
{
  if (inpoel.empty()) return;

  // find smallest node id
  auto minId = *std::min_element( begin(inpoel), end(inpoel) );

  // shift node ids to start from zero
  for (auto& n : inpoel) n -= minId;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsup( const std::vector< int >& inpoel, std::size_t nnpe )
//******************************************************************************
//  Generate derived data structure, elements surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< int > inpoel { 12, 14,  9, 11,
//!                                 10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \return Linked lists storing elements surrounding points
//! \warning It is not okay to call this function with an empty container or a
//!   non-positive number of nodes per element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _esup1_ and _esup2_, where _esup2_ holds the
//!   indices at which _esup1_ holds the element ids surrounding points. Looping
//!   over all elements surrounding all points can then be accomplished by the
//!   following loop:
//!   \code{.cpp}
//!     for (std::size_t p=0; p<npoin; ++p)
//!       for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i)
//!          use element id esup.first[i]
//!   \endcode
//!     To find out the number of points, _npoin_, the mesh connectivity,
//!     _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = static_cast< std::size_t >( *minmax.second + 1 );
//!   \endcode
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsup() with zero nodes per element" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = static_cast< std::size_t >( *minmax.second + 1 );

  // allocate one of the linked lists storing elements surrounding points: esup2
  // fill with zeros
  std::vector< std::size_t > esup2( npoin+1, 0 );

  // element pass 1: count number of elements connected to each point
  for (auto n : inpoel) ++esup2[ static_cast<std::size_t>(n) + 1 ];

  // storage/reshuffling pass 1: update storage counter and store
  // also find out the maximum size of esup1 (mesup)
  auto mesup = esup2[0]+1;
  for (std::size_t i=1; i<npoin+1; ++i) {
    esup2[i] += esup2[i-1];
    if (esup2[i]+1 > mesup) mesup = esup2[i]+1;
  }

  // now we know mesup, so allocate the other one of the linked lists storing
  // elements surrounding points: esup1
  std::vector< std::size_t > esup1( mesup );

  // store the elements in esup1
  std::size_t e = 0;
  for (auto n : inpoel) {
    auto N = static_cast< std::size_t >( n );
    auto j = esup2[N]+1;
    esup2[N] = j;
    esup1[j] = e/nnpe;
    ++e;
  }

  // storage/reshuffling pass 2
  auto snpoin = static_cast< long long >( npoin );
  for (auto i=snpoin; i>0; --i) {
    auto I = static_cast< std::size_t >( i );
    esup2[I] = esup2[I-1];
  }
  esup2[0] = 0;

  // Return (move out) linked list
  return std::make_pair( std::move(esup1), std::move(esup2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genPsup( const std::vector< int >& inpoel,
         std::size_t nnpe,
         const std::pair< std::vector< std::size_t >,
                          std::vector< std::size_t > >& esup )
//******************************************************************************
//  Generate derived data structure, points surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< int > inpoel { 12, 14,  9, 11,
//!                                 10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing points surrounding points
//! \warning It is not okay to call this function with an empty container or a
//!   non-positive number of nodes per element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _psup1_ and _psup2_, where _psup2_ holds the
//!   indices at which _psup1_ holds the point ids surrounding points. Looping
//!   over all points surrounding all points can then be accomplished by the
//!   following loop:
//!   \code{.cpp}
//!     for (std::size_t p=0; p<npoin; ++p)
//!       for (auto i=psup.second[p]+1; i<=psup.second[p+1]; ++i)
//!          use point id psup.first[i]
//!   \endcode
//!     To find out the number of points, _npoin_, the mesh connectivity,
//!     _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = static_cast< std::size_t >( *minmax.second + 1 );
//!   \endcode
//! \endcode
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genPsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genPsup() with zero nodes per element" );
  Assert( !esup.first.empty(), "Attempt to call genPsup() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genPsup() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = static_cast< std::size_t >( *minmax.second + 1 );

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // allocate both of the linked lists storing points surrounding points, we
  // only know the size of psup2, put in a single zero in psup1
  std::vector< std::size_t > psup2( npoin+1 );
  std::vector< std::size_t > psup1( 1, 0 );

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // fill both psup1 and psup2
  psup2[0] = 0;
  std::size_t j = 0;
  for (std::size_t p=0; p<npoin; ++p) {
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i ) {
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = static_cast< std::size_t >( inpoel[ esup1[i] * nnpe + n ] );
        if (q != p && lpoin[q] != p+1) {
          ++j;
          psup1.push_back( q );
          lpoin[q] = p+1;
        }
      }
    }
    psup2[p+1] = j;
  }

  // sort point ids for each point in psup1
  for (std::size_t p=0; p<npoin; ++p)
    std::sort(
      std::next( begin(psup1), static_cast<std::ptrdiff_t>(psup2[p]+1) ),
      std::next( begin(psup1), static_cast<std::ptrdiff_t>(psup2[p+1]+1) ) );

  // Return (move out) linked list
  return std::make_pair( std::move(psup1), std::move(psup2) );
}

} // tk::
