//******************************************************************************
/*!
  \file      src/Mesh/DerivedData.C
  \author    J. Bakosi
  \date      Fri 27 Mar 2015 01:46:55 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
//******************************************************************************

#include <algorithm>
#include <set>

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
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
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
  std::vector< std::size_t > psup2( npoin+1 ), psup1( 1, 0 );

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

std::vector< std::size_t >
genInpoed( const std::vector< int >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup )
//******************************************************************************
//  Generate derived data structure, edge connectivity
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
//! \return Linear vector storing node ids connecting edges
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linear vector with both
//!   point ids connecting an edge for all edges. By the name inpoed, Rainald
//!   probably meant interconnectivity of points and edges.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genInpoed() on empty container" );
  Assert( nnpe > 0, "Attempt to call genInpoed() with zero nodes per element" );
  Assert( !esup.first.empty(), "Attempt to call genInpoed() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genInpoed() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = static_cast< std::size_t >( *minmax.second + 1 );

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // map to contain stars, a point associated to points connected with edges
  std::map< std::size_t, std::vector< std::size_t > > star;

  // generate edge connectivity and store as stars where center id < spike id
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i )
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = static_cast< std::size_t >( inpoel[ esup1[i] * nnpe + n ] );
        if (q != p && lpoin[q] != p+1) {
          if (p < q) star[p].push_back(q);
          lpoin[q] = p+1;
        }
      }

  // linear vector to store edge connectivity
  std::vector< std::size_t > inpoed;

  // sort non-center points of each star and store all in linear vector
  for (auto p : star) {
    std::sort( begin(p.second), end(p.second) );
    for (const auto& e : p.second) {
      inpoed.push_back( p.first );
      inpoed.push_back( e );
    }
  }

  // Return (move out) edge connectivity
  return inpoed;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsupel( const std::vector< int >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup )
//******************************************************************************
//  Generate derived data structure, elements surrounding points of elements
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
//! \return Linked lists storing elements surrounding points of elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _esupel1_ and _esupel2_, where _esupel2_
//!   holds the indices at which _esupel1_ holds the element ids surrounding
//!   points of elements. Looping over all elements surrounding the points of
//!   all elements can then be accomplished by the following loop:
//!   \code{.cpp}
//!     for (std::size_t e=0; e<nelem; ++e)
//!       for (auto i=esupel.second[e]+1; i<=esupel.second[e+1]; ++i)
//!          use element id esupel.first[i]
//!   \endcode
//!     To find out the number of elements, _nelem_, the size of the mesh
//!     connectivity vector, _inpoel_, can be devided by the number of nodes per
//!     elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
//! \endcode
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsupel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsupel() with zero nodes per element" );
  Assert( !esup.first.empty(), "Attempt to call genEsupel() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genEsupel() with empty esup2" );

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // linked lists storing elements surrounding points of elements, put in a
  // single zero in both
  std::vector< std::size_t > esupel2( 1, 0 ), esupel1( 1, 0 );

  std::size_t e = 0;
  std::set< std::size_t > esuel;
  for (auto p : inpoel) {       // loop over all points of all elements
    auto P = static_cast< std::size_t >( p );
    // collect unique element ids of elements surrounding points of element
    for (auto i=esup2[P]+1; i<=esup2[P+1]; ++i) esuel.insert( esup1[i] );
    if (++e%nnpe == 0) {        // when finished checking all nodes of element
      // erase element whose surrounding elements are considered
      esuel.erase( e/nnpe-1 );
      // store unique element ids in esupel1
      for (auto i : esuel) esupel1.push_back( i );
      // store end-index for element used to address into esupel1
      esupel2.push_back( esupel2.back() + esuel.size() );
      esuel.clear();
    }
  }

  // Return (move out) linked list
  return std::make_pair( std::move(esupel1), std::move(esupel2) );
}

std::vector< long int >
genEsuel( const std::vector< int >& inpoel,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esupel )
//******************************************************************************
//  Generate derived data structure, elements surrounding elements (tets only)
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< int > inpoel { 12, 14,  9, 11,
//!                                 10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] esupel Elements surrounding elements as linked lists, see
//!   tk::genEsupel.
//! \return Linear vector storing all 4 elements surrounding elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esupel.first or esupel.second; it will throw an exception.
//! \details The data generated here is stored in a linear vector with all
//!   maximum 4 element ids surrounding elements. Elements surrounding elements
//!   are those that share a face with the element. This is similar to elements
//!   surrounding points of elements, computed by tk::genEsupel(), but esuel is
//!   smaller as it only stores the maximum 4 element ids that share a face with
//!   the given element. If a face does not have a neighbor element its value in
//!    esuel is -1.
//! \note This function only works for tetrahedra element connectivity.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsupel() on empty container" );
  Assert( inpoel.size()%4 == 0, "Size of inpoel must be divisible by four" );
  Assert( !esupel.first.empty(), "Attempt to call genEsupel() with empty esup1" );
  Assert( !esupel.second.empty(), "Attempt to call genEsupel() with empty esup2" );

  auto& esupel1 = esupel.first;
  auto& esupel2 = esupel.second;

  auto nelem = inpoel.size()/4;

  // linear vector storing elements surrounding elements, a priori known size,
  // initialize with -1 indicating no element on that side, i.e., boundary
  std::vector< long int > esuel( inpoel.size(), -1 );

  for (std::size_t e=0; e<nelem; ++e) {
    auto A = inpoel[e*4+0];
    auto B = inpoel[e*4+1];
    auto C = inpoel[e*4+2];
    auto D = inpoel[e*4+3];
    std::set< long int > faces; // will collect elem ids of shared faces
    // loop over elements surrounding points of elements
    for (auto i=esupel2[e]+1; i<=esupel2[e+1]; ++i) {
      auto f = esupel1[i];
      auto F = static_cast< long int >( f );
      // test if a face is shared between elements e and f
      if ( (A==inpoel[f*4+0] || A==inpoel[f*4+1] ||
            A==inpoel[f*4+2] || A==inpoel[f*4+3]) &&
           (B==inpoel[f*4+0] || B==inpoel[f*4+1] ||
            B==inpoel[f*4+2] || B==inpoel[f*4+3]) &&
           (C==inpoel[f*4+0] || C==inpoel[f*4+1] ||
            C==inpoel[f*4+2] || C==inpoel[f*4+3]) ) faces.insert( F );

      if ( (A==inpoel[f*4+0] || A==inpoel[f*4+1] ||
            A==inpoel[f*4+2] || A==inpoel[f*4+3]) &&
           (B==inpoel[f*4+0] || B==inpoel[f*4+1] ||
            B==inpoel[f*4+2] || B==inpoel[f*4+3]) &&
           (D==inpoel[f*4+0] || D==inpoel[f*4+1] ||
            D==inpoel[f*4+2] || D==inpoel[f*4+3]) ) faces.insert( F );

      if ( (A==inpoel[f*4+0] || A==inpoel[f*4+1] ||
            A==inpoel[f*4+2] || A==inpoel[f*4+3]) &&
           (C==inpoel[f*4+0] || C==inpoel[f*4+1] ||
            C==inpoel[f*4+2] || C==inpoel[f*4+3]) &&
           (D==inpoel[f*4+0] || D==inpoel[f*4+1] ||
            D==inpoel[f*4+2] || D==inpoel[f*4+3]) ) faces.insert( F );

      if ( (C==inpoel[f*4+0] || C==inpoel[f*4+1] ||
            C==inpoel[f*4+2] || C==inpoel[f*4+3]) &&
           (B==inpoel[f*4+0] || B==inpoel[f*4+1] ||
            B==inpoel[f*4+2] || B==inpoel[f*4+3]) &&
           (D==inpoel[f*4+0] || D==inpoel[f*4+1] ||
            D==inpoel[f*4+2] || D==inpoel[f*4+3]) ) faces.insert( F );

      // store element ids for shared faces
      std::size_t k=0;
      for (auto j : faces) esuel[e*4+k++] = j;
    }
  }

  // Return (move out) linked list
  return esuel;
}

} // tk::
