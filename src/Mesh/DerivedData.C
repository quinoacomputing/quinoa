//******************************************************************************
/*!
  \file      src/Mesh/DerivedData.C
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 06:34:21 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
//******************************************************************************

#include <cstring>

#include <make_unique.h>
#include <DerivedData.h>

namespace tk {

std::pair< std::unique_ptr< std::size_t[] >, std::unique_ptr< std::size_t[] > >
genEsup( const std::vector< std::size_t >& inpoel, std::size_t npoin )
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
  auto esup2 = tk::make_unique< std::size_t[] >( npoin+1 );

  // element pass 1: count number of elements connected to each point
  std::memset( esup2.get(), 0, (npoin+1)*sizeof(std::size_t) );
  for (const auto& e : inpoel) ++esup2[e+1];

  // storage/reshuffling pass 1: update storage counter and store
  // also find out the maximum size of esup1 (mesup)
  std::size_t mesup = esup2[0]+1;
  for (std::size_t i=1; i<npoin+1; ++i) {
    esup2[i] += esup2[i-1];
    if (esup2[i]+1 > mesup) mesup = esup2[i]+1;
  }

  // now we know mesup, so allocate the other one of the linked lists storing
  // elements surrounding points: esup1
  auto esup1 = tk::make_unique< std::size_t[] >( mesup );

  // store the elements in esup1
  std::size_t t=0;
  for (const auto& e : inpoel) {
    std::size_t j = esup2[e]+1;
    esup2[e] = j;
    esup1[j] = t;
    ++t;
  }

  // storage/reshuffling pass 2
  for (std::size_t i=npoin; i>0; --i) esup2[i] = esup2[i-1];
  esup2[0] = 0;

  // Return (move out) linked list
  return { std::move(esup1), std::move(esup2) };
}

} // tk::
