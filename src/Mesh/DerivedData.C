//******************************************************************************
/*!
  \file      src/Mesh/DerivedData.C
  \author    J. Bakosi
  \date      Tue 17 Mar 2015 07:49:17 AM MDT
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
//! \author J. Bakosi
//******************************************************************************
{
  // find smallest node id
  auto minId = *std::min_element( begin(inpoel), end(inpoel) );

  // shift node ids to start from zero
  for (auto& n : inpoel) n -= minId;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsup( const std::vector< int >& inpoel, std::size_t npoin )
//******************************************************************************
//  Generate derived data structure, elements surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements
//! \param[in] npoin Number of points in mesh
//! \return Linked lists storing elements surrounding points
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  // allocate one of the linked lists storing elements surrounding points: esup2
  // fill with zeros
  std::vector< std::size_t > esup2( npoin+1, 0 );

  // element pass 1: count number of elements connected to each point
  for (auto e : inpoel) ++esup2[ static_cast<std::size_t>(e)+1 ];

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
  std::size_t t=0;
  for (auto e : inpoel) {
    auto E = static_cast< std::size_t >( e );
    auto j = esup2[E]+1;
    esup2[E] = j;
    esup1[j] = t;
    ++t;
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

} // tk::
