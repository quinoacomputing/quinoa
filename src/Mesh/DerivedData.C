// *****************************************************************************
/*!
  \file      src/Mesh/DerivedData.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
// *****************************************************************************

#include <set>
#include <map>
#include <iterator>
#include <algorithm>
#include <type_traits>
#include <cstddef>
#include <array>
#include <unordered_set>

#include "Exception.h"
#include "DerivedData.h"
#include "ContainerUtil.h"
#include "Vector.h"

namespace tk {

std::size_t
npoin( const std::vector< std::size_t >& inpoel )
// *****************************************************************************
// Compute number of points (nodes) in mesh from connectivity
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//! \return Number of mesh points (nodes)
// *****************************************************************************
{
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  return *minmax.second + 1;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsup( const std::vector< std::size_t >& inpoel, std::size_t nnpe )
// *****************************************************************************
//  Generate derived data structure, elements surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
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
//!     auto npoin = *minmax.second + 1;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsup() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // allocate one of the linked lists storing elements surrounding points: esup2
  // fill with zeros
  std::vector< std::size_t > esup2( npoin+1, 0 );

  // element pass 1: count number of elements connected to each point
  for (auto n : inpoel) ++esup2[ n + 1 ];

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
    auto j = esup2[n]+1;
    esup2[n] = j;
    esup1[j] = e/nnpe;
    ++e;
  }

  // storage/reshuffling pass 2
  for (auto i=npoin; i>0; --i) esup2[i] = esup2[i-1];
  esup2[0] = 0;

  // Return (move out) linked lists
  return std::make_pair( std::move(esup1), std::move(esup2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genPsup( const std::vector< std::size_t >& inpoel,
         std::size_t nnpe,
         const std::pair< std::vector< std::size_t >,
                          std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, points surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
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
//!    To find out the number of points, _npoin_, the mesh connectivity,
//!    _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = *minmax.second + 1;
//!   \endcode
//!   or the length-1 of the generated index list:
//!   \code{.cpp}
//!     auto npoin = psup.second.size()-1;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genPsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genPsup() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genPsup() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genPsup() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

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
        auto q = inpoel[ esup1[i] * nnpe + n ];
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

  // Return (move out) linked lists
  return std::make_pair( std::move(psup1), std::move(psup2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEdsup( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, edges surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element (3 or 4)
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing edges (point ids p < q) emanating from points
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _edsup1_ and _edsup2_, where _edsup2_ holds
//!   the indices at which _edsup1_ holds the edge-end point ids emanating from
//!   points for all points. The generated data structure, linked lists edsup1
//!   and edsup2, are very similar to psup1 and psup2, generated by genPsup(),
//!   except here only unique edges are stored, i.e., for edges with point ids
//!   p < q, only ids q are stored that are still associated to point p. Looping
//!   over all unique edges can then be accomplished by the following loop:
//!   \code{.cpp}
//!     for (std::size_t p=0; p<npoin; ++p)
//!       for (auto i=edsup.second[p]+1; i<=edsup.second[p+1]; ++i)
//!         use edge with point ids p < edsup.first[i]
//!   \endcode
//!   To find out the number of points, _npoin_, the mesh connectivity,
//!   _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = *minmax.second + 1;
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see tk::genInpoed for similar data that sometimes may be more advantageous
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEdsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEdsup() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genEdsup() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genEdsup() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genEdsup() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // map to contain stars, a point associated to points connected with edges
  // storing only the end-point id, q, of point ids p < q
  std::map< std::size_t, std::vector< std::size_t > > star;

  // generate edge connectivity and store as stars where center id < spike id
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i )
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = inpoel[ esup1[i] * nnpe + n ];
        if (q != p && lpoin[q] != p+1) {
          if (p < q) star[p].push_back(q);
          lpoin[q] = p+1;
        }
      }

  // linked lists (vectors) to store edges surrounding points and their indices
  std::vector< std::size_t > edsup1( 1, 0 ), edsup2( 1, 0 );

  // sort non-center points of each star and store nodes and indices in vectors
  for (auto p : star) {
    std::sort( begin(p.second), end(p.second) );
    edsup2.push_back( edsup2.back() + p.second.size() );
    for (auto e : p.second) edsup1.push_back( e );
  }
  // fill up index array with the last index for points with no new edges
  for (std::size_t i=0; i<npoin-star.size(); ++i)
    edsup2.push_back( edsup2.back() );

  // Return (move out) linked lists
  return std::make_pair( std::move(edsup1), std::move(edsup2) );
}

std::vector< std::size_t >
genInpoed( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, edge connectivity
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element (3 or 4)
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linear vector storing edge connectivity (point ids p < q)
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linear vector and is very
//!   similar to the linked lists, _edsup1_ and _edsup2, generated by
//!   genEdsup(). The difference is that in the linear vector, inpoed, generated
//!   here, both edge point ids are stored as a pair, p < q, as opposed to the
//!   linked lists edsup1 and edsup2, in which edsup1 only stores the edge-end
//!   point ids (still associated to edge-start point ids when used together
//!   with edsup2). The rationale is that while inpoed is larger in memory, it
//!   allows direct access to edges (pair of point ids making up an edge),
//!   edsup1 and edsup2 are smaller in memory, still allow accessing the same
//!   data (edge point id pairs) but only in a linear fashion, not by direct
//!   access to particular edges. Accessing all unique edges using the edge
//!   connectivity data structure, inpoed, generated here can be accomplished by
//!   \code{.cpp}
//!     for (std::size_t e=0; e<inpoed.size()/2; ++e) {
//!       use point id p of edge e = inpoed[e*2];
//!       use point id q of edge e = inpoed[e*2+1];
//!     }
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see tk::genEdsup for similar data that sometimes may be more advantageous
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genInpoed() on empty container" );
  Assert( nnpe > 0, "Attempt to call genInpoed() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genInpoed() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genInpoed() with empty esup1" );
  Assert( !esup.second.empty(),
          "Attempt to call genInpoed() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // map to contain stars, a point associated to points connected with edges,
  // storing only the end-point id, q, of point ids p < q
  std::map< std::size_t, std::vector< std::size_t > > star;

  // generate edge connectivity and store as stars where center id < spike id
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i )
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = inpoel[ esup1[i] * nnpe + n ];
        if (q != p && lpoin[q] != p+1) {
          if (p < q) star[p].push_back( q );
          lpoin[q] = p+1;
        }
      }

  // linear vector to store edge connectivity and their indices
  std::vector< std::size_t > inpoed;

  // sort non-center points of each star and store both start and end points of
  // each star in linear vector
  for (auto p : star) {
    std::sort( begin(p.second), end(p.second) );
    for (auto e : p.second) {
      inpoed.push_back( p.first );
      inpoed.push_back( e );
    }
  }

  // Return (move out) linear vector
  return inpoed;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsupel( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding points of elements
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
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
//!   To find out the number of elements, _nelem_, the size of the mesh
//!   connectivity vector, _inpoel_, can be devided by the number of nodes per
//!   elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsupel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsupel() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genEsupel() with empty esup1" );
  Assert( !esup.second.empty(),
          "Attempt to call genEsupel() with empty esup2" );

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // linked lists storing elements surrounding points of elements, put in a
  // single zero in both
  std::vector< std::size_t > esupel2( 1, 0 ), esupel1( 1, 0 );

  std::size_t e = 0;
  std::set< std::size_t > esuel;
  for (auto p : inpoel) {       // loop over all points of all elements
    // collect unique element ids of elements surrounding points of element
    for (auto i=esup2[p]+1; i<=esup2[p+1]; ++i) esuel.insert( esup1[i] );
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

  // Return (move out) linked lists
  return std::make_pair( std::move(esupel1), std::move(esupel2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsuel( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding elements
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing elements surrounding elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _esuel1_ and _esuel2_, where _esuel2_ holds
//!   the indices at which _esuel1_ holds the element ids surrounding elements.
//!   Looping over elements surrounding elements can then be accomplished by the
//!   following loop:
//!   \code{.cpp}
//!     for (std::size_t e=0; e<nelem; ++e)
//!       for (auto i=esuel.second[e]+1; i<=esuel.second[e+1]; ++i)
//!          use element id esuel.first[i]
//!   \endcode
//!   To find out the number of elements, _nelem_, the size of the mesh
//!   connectivity vector, _inpoel_, can be devided by the number of nodes per
//!   elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
//! \note In principle, this function *should* work for any positive nnpe,
//!   however, only nnpe = 4 (tetrahedra) and nnpe = 3 (triangles) are tested.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsuel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsuel() with zero nodes per element" );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by four" );
  Assert( !esup.first.empty(), "Attempt to call genEsuel() with empty esuel1" );
  Assert( !esup.second.empty(),
          "Attempt to call genEsuel() with empty esuel2" );

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  auto nelem = inpoel.size()/nnpe;

  // lambda that returns true if elements hel and gel share a face
  auto adj = [ &inpoel, nnpe ]( std::size_t hel, std::size_t gel ) -> bool {
    std::vector< bool > sp;
    for (std::size_t h=0; h<nnpe; ++h)
      for (std::size_t g=0; g<nnpe; ++g)
        if (inpoel[hel*nnpe+h] == inpoel[gel*nnpe+g]) sp.push_back( true );
    if (sp.size() == nnpe-1) return true; else return false;
  };

  // map to associate unique elements and their surrounding elements
  std::map< std::size_t, std::vector< std::size_t > > es;

  for (std::size_t e=0; e<nelem; ++e) {
    std::set< std::size_t > faces; // will collect elem ids of shared faces
    for (std::size_t n=0; n<nnpe; ++n) {
      auto i = inpoel[ e*nnpe+n ];
      for (auto j=esup2[i]+1; j<=esup2[i+1]; ++j)
        if (adj( e, esup1[j] )) faces.insert( esup1[j] );
    }
    // store element ids of shared faces
    for (auto j : faces) es[e].push_back(j);
  }

  // storing elements surrounding elements
  std::vector< std::size_t > esuel1( 1, 0 ), esuel2( 1, 0 );

  // store elements surrounding elements in linked lists
  for (auto e : es) {
    esuel2.push_back( esuel2.back() + e.second.size() );
    for (auto s : e.second) esuel1.push_back( s );
  }

  // Return (move out) linked lists
  return std::make_pair( std::move(esuel1), std::move(esuel2) );
}

std::vector< std::size_t >
genInedel( const std::vector< std::size_t >& inpoel,
           std::size_t nnpe,
           const std::vector< std::size_t >& inpoed )
// *****************************************************************************
//  Generate derived data structure, edges of elements
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element
//! \param[in] inpoed Edge connectivity as linear vector, see tk::genInpoed
//! \return Linear vector storing all edge ids * 2 of all elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or inpoed or a non-positive number of nodes per element; it will
//!   throw an exception.
//! \details The data generated here is stored in a linear vector with all
//!   edge ids (as defined by inpoed) of all elements. The edge ids stored in
//!   inedel can be directly used to index the vector inpoed. Because the
//!   derived data structure generated here, inedel, is intended to be used in
//!   conjunction with the linear vector inpoed and not with the linked lists
//!   edsup1 and edsup2, this function takes inpoed as an argument. Accessing
//!   the edges of element e using the edge of elements data structure, inedel,
//!   generated here can be accomplished by
//!   \code{.cpp}
//!     for (std::size_t e=0; e<nelem; ++e) {
//!       for (std::size_t i=0; i<nepe; ++i) {
//!         use edge id inedel[e*nepe+i] of element e, or
//!         use point ids p < q of edge id inedel[e*nepe+i] of element e as
//!           p = inpoed[ inedel[e*nepe+i]*2 ]
//!           q = inpoed[ inedel[e*nepe+i]*2+1 ]
//!       }
//!     }
//!   \endcode
//!   where _nepe_ denotes the number of edges per elements: 3 for triangles, 6
//!   for tetrahedra. To find out the number of elements, _nelem_, the size of
//!   the mesh connectivity vector, _inpoel_, can be devided by the number of
//!   nodes per elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genInedel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genInedel() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genInedel() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !inpoed.empty(), "Attempt to call genInedel() with empty inpoed" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // First, generate index of star centers. This is necessary to avoid a
  // brute-force search for point ids of edges when searching for element edges.
  // Note that this is the same as edsup2, generated by genEdsup(). However,
  // because the derived data structure generated here, inedel, is intended to
  // be used in conjunction with the linear vector inpoed and not with the
  // linked lists edsup1 and edsup2, this function takes inpoed as an argument,
  // and so edsup2 is temporarily generated here to avoid a brute-force search.

  // map to contain stars, a point associated to points connected with edges
  // storing only the end-point id, q, of point ids p < q
  std::map< std::size_t, std::vector< std::size_t > > star;

  // generate stars from inpoed; starting with zero, every even is a star
  // center, every odd is a spike
  for (std::size_t i=0; i<inpoed.size()/2; ++i)
    star[ inpoed[i*2] ].push_back( inpoed[i*2+1] );

  // store index of star centers in vector; assume non-center points of each
  // star have already been sorted
  std::vector< std::size_t > edsup2( 1, 0 );
  for (auto p : star) edsup2.push_back( edsup2.back() + p.second.size() );
  // fill up index array with the last index for points with no new edges
  for (std::size_t i=0; i<npoin-star.size(); ++i)
    edsup2.push_back( edsup2.back() );
  star.clear();

  // Second, generate edges of elements

  auto nelem = inpoel.size()/nnpe;

  // map associating elem id with vector of edge ids
  std::map< std::size_t, std::vector< std::size_t > > edges;

  // generate map of elements associated to edge ids
  for (std::size_t e=0; e<nelem; ++e)
    for (std::size_t n=0; n<nnpe; ++n) {
      auto p = inpoel[e*nnpe+n];
      for (auto i=edsup2[p]+1; i<=edsup2[p+1]; ++i)
         for (std::size_t j=0; j<nnpe; ++j)
            if (inpoed[(i-1)*2+1] == inpoel[e*nnpe+j])
              edges[e].push_back( i-1 );
    }

  // linear vector to store the edge ids of all elements
  std::vector< std::size_t > inedel;

  // store edge ids of elements in linear vector
  for (auto e : edges) for (auto p : e.second) inedel.push_back( p );

  // Return (move out) vector
  return inedel;
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEsued( const std::vector< std::size_t >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding edges
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] nnpe Number of nodes per element (3 or 4)
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Linked lists storing elements surrounding edges
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second or a non-positive number of nodes per
//!   element; it will throw an exception.
//! \details The data generated here is stored in a linked list, more precisely,
//!   two linked arrays (vectors), _esued1_ and _esued2_, where _esued2_ holds
//!   the indices at which _esued1_ holds the element ids surrounding edges.
//!   Looping over all elements surrounding edges can then be accomplished by
//!   the following loop:
//!   \code{.cpp}
//!     for (std::size_t e=0; e<nedge; ++e)
//!       for (auto i=esued.second[e]+1; i<=esued.second[e+1]; ++i)
//!         use element id esued.first[i]
//!   \endcode
//!   To find out the number of edges, _nedge_, the edge connectivity, _inpoed_,
//!   can be queried:
//!   \code{.cpp}
//!     auto esup = tk::genEsup(inpoel,nnpe);
//!     auto nedge = tk::genInpoed(inpoel,nnpe,esup).size()/2;
//!   \endcode
//!   where _nnpe_ is the number of nodes per element (4 for tetrahedra, 3 for
//!   triangles).
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsued() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEsued() with zero nodes per element" );
  Assert( nnpe == 3 || nnpe == 4,
          "Attempt to call genEsued() with nodes per element, nnpe, that is "
          "neither 4 (tetrahedra) nor 3 (triangles)." );
  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by nnpe" );
  Assert( !esup.first.empty(), "Attempt to call genEsued() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genEsued() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // allocate and fill with zeros a temporary array, only used locally
  std::vector< std::size_t > lpoin( npoin, 0 );

  // lambda that returns true if element e contains edge (p < q)
  auto has = [ &inpoel, nnpe ]( std::size_t e, std::size_t p, std::size_t q )
  -> bool {
    std::vector< bool > sp;
    for (std::size_t n=0; n<nnpe; ++n)
      if (inpoel[e*nnpe+n] == p || inpoel[e*nnpe+n] == q)
        sp.push_back( true );
    if (sp.size() == 2) return true; else return false;
  };

  // map to associate edges to unique surrounding element ids
  std::map< std::size_t,  std::vector< std::size_t > > revolver;

  // generate edges and associated vector of unique surrounding element ids
  std::size_t ed = 0;
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t i=esup2[p]+1; i<=esup2[p+1]; ++i )
      for (std::size_t n=0; n<nnpe; ++n) {
        auto q = inpoel[ esup1[i] * nnpe + n ];
        if (q != p && lpoin[q] != p+1) {
          if (p < q) {  // for edge given point ids p < q
            for (std::size_t j=esup2[p]+1; j<=esup2[p+1]; ++j ) {
              auto e = esup1[j];
              if (has(e,p,q)) revolver[ed].push_back(e);
            }
            ++ed;
          }
          lpoin[q] = p+1;
        }
      }

  // linked lists (vectors) to store elements surrounding edges
  std::vector< std::size_t > esued1( 1, 0 ), esued2( 1, 0 );

  // sort and store elements surrounding edges and their indices in vectors
  for (auto p : revolver) {
    std::sort( begin(p.second), end(p.second) );
    esued2.push_back( esued2.back() + p.second.size() );
    for (auto e : p.second) esued1.push_back( e );
  }

  // Return (move out) linked lists
  return std::make_pair( std::move(esued1), std::move(esued2) );
}

std::size_t
genNbfacTet( std::size_t tnbfac,
             const std::vector< std::size_t >& inpoel,
             const std::vector< std::size_t >& triinpoel_complete,
             const std::map< int, std::vector< std::size_t > >& bface_complete,
             const std::unordered_map< std::size_t, std::size_t >& lid,
             std::vector< std::size_t >& triinpoel,
             std::map< int, std::vector< std::size_t > >& bface )
// *****************************************************************************
//  Generate the number of boundary-faces and the triangle boundary-face
//  connectivity for a chunk of a full mesh.
//  \warning This is for Triangular face-elements only.
//! \param[in] tnbfac Total number of boundary faces in the entire mesh.
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh.
//! \param[in] triinpoel_complete Interconnectivity of points and boundary-face
//!   in the entire mesh.
//! \param[in] bface_complete Map of boundary-face lists mapped to corresponding 
//!   side set ids for the entire mesh.
//! \param[in] lid Mapping between the node indices used in the smaller inpoel
//!   connectivity (a subset of the entire triinpoel_complete connectivity),
//!   e.g., after mesh partitioning.
//! \param[inout] triinpoel Interconnectivity of points and boundary-face in
//!   this mesh-partition.
//! \param[inout] bface Map of boundary-face lists mapped to corresponding 
//!   side set ids for this mesh-partition
//! \return Number of boundary-faces on this chare/mesh-partition.
//! \details This function takes a mesh by its domain-element
//!   (tetrahedron-connectivity) in inpoel and a boundary-face (triangle)
//!   connectivity in triinpoel_complete. Based on these two arrays, it
//!   searches for those faces of triinpoel_complete that are also in inpoel
//!   and as a result it generates (1) the number of boundary faces shared with
//!   the mesh in inpoel and (2) the intersection of the triangle element
//!   connectivity whose faces are shared with inpoel. An example use case is
//!   where triinpoel_complete contains the connectivity for the boundary of the
//!   full problem/mesh and inpoel contains the connectivity for only a chunk of
//!   an already partitioned mesh. This function then intersects
//!   triinpoel_complete with inpoel and returns only those faces that share
//!   nodes with inpoel.
// *****************************************************************************
{
  std::size_t nbfac(0), nnpf(3);

  if (tnbfac > 0)
  {

  Assert( !inpoel.empty(),
          "Attempt to call genNbfacTet() on empty inpoel container" );
  Assert( !triinpoel_complete.empty(),
          "Attempt to call genNbfacTet() on empty triinpoel_complete container" );
  Assert( triinpoel_complete.size()/nnpf == tnbfac, 
          "Incorrect size of triinpoel in genNbfacTet()" );

  auto nptet = inpoel;
  auto nptri = triinpoel_complete;

  tk::unique( nptet );
  tk::unique( nptri );

  std::unordered_set< std::size_t > snptet;

  // getting the reduced inpoel as a set for quick searches
  snptet.insert( begin(nptet), end(nptet));

  // vector to store boundary-face-nodes in this chunk
  std::vector< std::size_t > nptri_chunk;

  // getting the nodes of the boundary-faces in this chunk
  for (auto i : nptri)
    if (snptet.find(i) != end(snptet))
      nptri_chunk.push_back(i);

  std::size_t tag, icoun;

  // matching nodes in nptri_chunk with nodes in inpoel and 
  // triinpoel_complete to get the number of faces in this chunk
  for (const auto& ss : bface_complete)
  {
    for (auto f : ss.second)
    {
      icoun = f*nnpf;
      tag = 0;
      for (std::size_t i=0; i<nnpf; ++i)
      {
        for (auto j : nptri_chunk)
        {
          if (triinpoel_complete[icoun+i] == j) ++tag;
        }
      }
      if (tag == nnpf)
      // this is a boundary face
      {
        for (std::size_t i=0; i<nnpf; ++i)
        {
          auto ip = triinpoel_complete[icoun+i];

          // find local renumbered node-id to store in triinpoel
          triinpoel.push_back( tk::cref_find(lid,ip) );
        }

        bface[ss.first].push_back(nbfac);
        ++nbfac;
      }
    }
  }

  }

  return nbfac;
}

std::vector< int >
genEsuelTet( const std::vector< std::size_t >& inpoel,
             const std::pair< std::vector< std::size_t >,
                              std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data structure, elements surrounding elements
//  as a fixed length data structure as a full vector, including
//  boundary elements as -1.
//  \warning This is for Tetrahedra only.
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< std::size_t > inpoel { 12, 14,  9, 11,
//!                                         10, 14, 13, 12 };
//!   \endcode
//!   specifies two tetrahedra whose vertices (node ids) are { 12, 14, 9, 11 },
//!   and { 10, 14, 13, 12 }.
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Vector storing elements surrounding elements
//! \warning It is not okay to call this function with an empty container for
//!   inpoel or esup.first or esup.second; it will throw an exception.
//! \details The data generated here is stored in a single vector, with length
//!   nfpe * nelem. Note however, that nelem is not explicitly provided, but
//!   calculated from inpoel. For boundary elements, at the boundary face, this
//!   esuelTet stores value -1 indicating that this is outside the domain. The
//!   convention for numbering the local face (triangle) connectivity is very
//!   important, e.g., in generating the inpofa array later. This node ordering
//!   convention is stored in tk::lpofa. Thus function is specific to
//!   tetrahedra, which is reflected in the fact that nnpe and nfpe are being
//!   set here in the function rather than being input arguments. To find out
//!   the number of elements, _nelem_, the size of the mesh connectivity vector,
//!   _inpoel_, can be devided by the number of nodes per elements, _nnpe_:
//!   \code{.cpp}
//!     auto nelem = inpoel.size()/nnpe;
//!   \endcode
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEsuelTet() on empty container" );
  Assert( !esup.first.empty(), "Attempt to call genEsuelTet() with empty esup1" );
  Assert( !esup.second.empty(),
          "Attempt to call genEsuelTet() with empty esup2" );

  auto& esup1 = esup.first;
  auto& esup2 = esup.second;

  // set tetrahedron geometry
  std::size_t nnpe(4), nfpe(4), nnpf(3);

  Assert( inpoel.size()%nnpe == 0, "Size of inpoel must be divisible by four" );

  // get nelem and npoin
  auto nelem = inpoel.size()/nnpe;
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  std::vector< int > esuelTet(nfpe*nelem, -1);
  std::vector< std::size_t > lhelp(nnpf,0),
                             lpoin(npoin,0);

  for (std::size_t e=0; e<nelem; ++e)
  {
    auto mark = nnpe*e;
    for (std::size_t fe=0; fe<nfpe; ++fe)
    {
      // array which stores points on this face
      lhelp[0] = inpoel[mark+lpofa[fe][0]];
      lhelp[1] = inpoel[mark+lpofa[fe][1]];
      lhelp[2] = inpoel[mark+lpofa[fe][2]];

      // mark in this array
      lpoin[lhelp[0]] = 1;
      lpoin[lhelp[1]] = 1;
      lpoin[lhelp[2]] = 1;

      // select a point on this face
      auto ipoin = lhelp[0];

      // loop over elements around this point
      for (std::size_t j=esup2[ipoin]+1; j<=esup2[ipoin+1]; ++j )
      {
        auto jelem = esup1[j];
        // if this jelem is not e itself then proceed
        if (jelem != e)
        {
          for (std::size_t fj=0; fj<nfpe; ++fj)
          {
            std::size_t icoun(0);
            for (std::size_t jnofa=0; jnofa<nnpf; ++jnofa)
            {
              auto markj = jelem*nnpe;
              auto jpoin = inpoel[markj+lpofa[fj][jnofa]];
              if (lpoin[jpoin] == 1) { ++icoun; }
            }
            //store esuel if
            if (icoun == nnpf)
            {
              auto markf = nfpe*e;
              esuelTet[markf+fe] = static_cast<int>(jelem);

              markf = nfpe*jelem;
              esuelTet[markf+fj] = static_cast<int>(e);
            }
          }
        }
      }
      // reset this array
      lpoin[lhelp[0]] = 0;
      lpoin[lhelp[1]] = 0;
      lpoin[lhelp[2]] = 0;
    }
  }

  return esuelTet;
}

std::size_t
genNtfac( std::size_t nfpe,
          std::size_t nbfac,
          const std::vector< int >& esuelTet )
// *****************************************************************************
//  Generate derived data, total number of faces in the mesh
//! \param[in] nfpe Number of faces per element.
//! \param[in] nbfac Number of boundary faces.
//! \param[in] esuelTet Elements surrounding elements.
//! \return Total number of faces in the mesh
//! \details The unsigned integer here gives the total number of faces in
//     the mesh.
// *****************************************************************************
{
  Assert( !esuelTet.empty(), "Attempt to call genNtfac() with empty esuelTet" );
  Assert( esuelTet.size()%nfpe == 0,
                  "Size of esuelTet must be divisible by nfpe" );
  Assert( nfpe > 0, "Attempt to call genNtfac() with zero faces per element" );

  auto nelem = esuelTet.size()/nfpe;

  std::size_t nifac = 0;

  if (nbfac > 0) {
    // loop through elements surrounding elements to find number of internal faces
    for (std::size_t e=0; e<nelem; ++e)
    {
      for (std::size_t ip=nfpe*e; ip<nfpe*(e+1); ++ip)
      {
        if (esuelTet[ip] != -1)
        {
          if ( e<static_cast< std::size_t >(esuelTet[ip]) )
          {
            ++nifac;
          }
        }
      }
    }
  }

  return nifac + nbfac;
}

std::vector< int >
genEsuf( std::size_t nfpe,
         std::size_t ntfac,
         std::size_t nbfac,
         const std::vector< std::size_t >& belem,
         const std::vector< int >& esuelTet )
// *****************************************************************************
//  Generate derived data, elements surrounding faces
//! \param[in] nfpe  Number of faces per element.
//! \param[in] ntfac Total number of faces.
//! \param[in] nbfac Number of boundary faces.
//! \param[in] belem Boundary element vector.
//! \param[in] esuelTet Elements surrounding elements.
//! \return Elements surrounding faces.
//! \details The unsigned integer vector gives the IDs of the elements to the
//    left and the right of each face in the mesh. The convention followed 
//    throughout is : The left element always has an ID smaller than the ID of
//    the right element.
// *****************************************************************************
{
  Assert( esuelTet.size()%nfpe == 0, 
                  "Size of esuelTet must be divisible by nfpe" );
  Assert( nfpe > 0, "Attempt to call genEsuf() with zero faces per element" );

  auto nelem = esuelTet.size()/nfpe;

  std::vector< int > esuf(2*ntfac);

  if (nbfac > 0) {
    // counters for number of internal and boundary faces
    std::size_t icoun(2*nbfac), bcoun(0);

    // loop to get face-element connectivity for internal faces
    for (std::size_t e=0; e<nelem; ++e) {
      for (std::size_t ip=nfpe*e; ip<nfpe*(e+1); ++ip) {
        auto jelem = esuelTet[ip];
        if (jelem != -1)
        {
          if ( e < static_cast< std::size_t >(jelem) )
          {
            esuf[icoun] = static_cast< int >(e);
            esuf[icoun+1] = static_cast< int >(jelem);
            icoun = icoun + 2;
          }
        }
      }
    }

    bcoun = 0;
    for (auto ie : belem) {
      esuf[bcoun] = static_cast< int >(ie);
      esuf[bcoun+1] = -1;  // outside domain
      bcoun = bcoun + 2;
    }
  }

  return esuf;
}

std::vector< std::size_t >
genInpofaTet( std::size_t ntfac,
              std::size_t nbfac,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::size_t >& triinpoel,
              const std::vector< int >& esuelTet )
// *****************************************************************************
//  Generate derived data, points on faces for tetrahedra only
//! \param[in] ntfac Total number of faces.
//! \param[in] nbfac Number of boundary faces.
//! \param[in] inpoel Element-node connectivity.
//! \param[in] triinpoel Face-node connectivity.
//! \param[in] esuelTet Elements surrounding elements.
//! \return Points surrounding faces. The unsigned integer vector gives the
//!   elements to the left and to the right of each face in the mesh.
// *****************************************************************************
{
  std::vector< std::size_t > inpofa;

  if (nbfac > 0)
  {
    // set tetrahedron geometry
    std::size_t nnpe(4), nfpe(4), nnpf(3);

    Assert( esuelTet.size()%nfpe == 0,
                    "Size of esuelTet must be divisible by nfpe" );
    Assert( inpoel.size()%nnpe == 0,
                    "Size of inpoel must be divisible by nnpe" );

    auto nelem = inpoel.size()/nnpe;

    inpofa.resize(nnpf*ntfac);

    // counters for number of internal and boundary faces
    std::size_t icoun(nnpf*nbfac);
    std::size_t mark(0);

    // loop over elems to get nodes on faces
    // this fills the interior face-node connectivity part
    for (std::size_t e=0; e<nelem; ++e)
    {
      mark = nnpe*e;
      for (std::size_t f=0; f<nfpe ; ++f)
      {
        auto ip = nfpe*e + f;
        auto jelem = esuelTet[ip];
        if (jelem != -1)
        {
          if ( e < static_cast< std::size_t >(jelem) )
          {
            inpofa[icoun]   = inpoel[mark+lpofa[f][0]];
            inpofa[icoun+1] = inpoel[mark+lpofa[f][1]];
            inpofa[icoun+2] = inpoel[mark+lpofa[f][2]];
            icoun = icoun + nnpf;
          }
        }
      }
    }

    // this fills the boundary face-node connectivity part
    // consistent with triinpoel
    for (std::size_t f=0; f<nbfac; ++f)
    {
      icoun = nnpf * f;
      inpofa[icoun+0] = triinpoel[icoun+2];
      inpofa[icoun+1] = triinpoel[icoun+1];
      inpofa[icoun+2] = triinpoel[icoun+0];
    }
  }

  return inpofa;
}
        
std::vector< std::size_t >
genBelemTet( std::size_t nbfac,
              const std::vector< std::size_t >& inpofa,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& esup )
// *****************************************************************************
//  Generate derived data, and array of elements which share one or more of
//   their faces with the domain boundary (host elements).
//! \param[in] nbfac Number of boundary faces.
//! \param[in] inpofa Face-node connectivity.
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Host elements or boundary elements. The unsigned integer vector
//!   gives the elements to the left of each boundary face in the mesh.
// *****************************************************************************
{
  std::vector< std::size_t > belem(nbfac);

  if (nbfac > 0)
  {

  // set tetrahedron geometry
  std::size_t nnpf(3), tag(0);

  // loop over all the boundary faces
  for(std::size_t f=0; f<nbfac; ++f)
  {
    belem[f] = 0;

    // array storing the element-cluster around face
    std::vector< std::size_t > elemcluster;

    // loop over the nodes of this boundary face
    for(std::size_t lp=0; lp<nnpf; ++lp)
    {
      auto gp = inpofa[nnpf*f + lp];

      Assert( gp < esup.second.size(), "Indexing out of esup2" );
      // loop over elements surrounding this node
      for (auto i=esup.second[gp]+1; i<=esup.second[gp+1]; ++i)
      {
        // form element-cluster vector
        elemcluster.push_back(esup.first[i]);
      }
    }

    // loop over element cluster to find repeating elements
    for(std::size_t i=0; i<elemcluster.size(); ++i)
    {
      auto ge = elemcluster[i];
      tag = 1;
      for(std::size_t j=0; j<elemcluster.size(); ++j)
      {
        if ( i != j && elemcluster[j] == ge )
        {
          tag++;
        }
      }
      if (tag == nnpf)
      {
        // this is the required boundary element
        belem[f] = ge;
        break;
      }
    }
  }
  }

  return belem;
}
        
tk::Fields
genGeoFaceTri( std::size_t ntfac,
               const std::vector< std::size_t >& inpofa,
               const tk::UnsMesh::Coords& coord )
// *****************************************************************************
//  Generate derived data, which stores the geometry details both internal and
//   boundary triangular faces in the mesh.
//! \param[in] ntfac Total number of faces in the mesh.
//! \param[in] inpofa Face-node connectivity.
//! \param[in] coord Co-ordinates of nodes in this mesh-chunk.
//! \return Face geometry information. This includes face area, unit normal
//!   pointing outward of the element to the left of the face, and face
//!   centroid coordinates. Use the following examples to access this
//!   information for face-f.
//!   face area: geoFace(f,0,0),
//!   unit-normal x-component: geoFace(f,1,0),
//!               y-component: geoFace(f,2,0),
//!               z-component: geoFace(f,3,0),
//!   centroid x-coordinate: geoFace(f,4,0),
//!            y-coordinate: geoFace(f,5,0),
//!            z-coordinate: geoFace(f,6,0).
// *****************************************************************************
{
  tk::Fields geoFace( ntfac, 7 );

  // set triangle geometry
  std::size_t nnpf(3);

  Assert( inpofa.size()%nnpf == 0,
          "Size of inpofa must be divisible by nnpf" );

  for(std::size_t f=0; f<ntfac; ++f)
  {
    std::size_t ip1, ip2, ip3;
    tk::real xp1, yp1, zp1,
             xp2, yp2, zp2,
             xp3, yp3, zp3;

    // get area
    ip1 = inpofa[nnpf*f];
    ip2 = inpofa[nnpf*f + 1];
    ip3 = inpofa[nnpf*f + 2];

    xp1 = coord[0][ip1];
    yp1 = coord[1][ip1];
    zp1 = coord[2][ip1];

    xp2 = coord[0][ip2];
    yp2 = coord[1][ip2];
    zp2 = coord[2][ip2];

    xp3 = coord[0][ip3];
    yp3 = coord[1][ip3];
    zp3 = coord[2][ip3];

    auto geoif = geoFaceTri( {{xp1, xp2, xp3}},
                             {{yp1, yp2, yp3}},
                             {{zp1, zp2, zp3}} );

    for (std::size_t i=0; i<7; ++i)
      geoFace(f,i,0) = geoif(0,i,0);
  }

  return geoFace;
}

tk::Fields
geoFaceTri( const std::array< tk::real, 3 >& x,
            const std::array< tk::real, 3 >& y,
            const std::array< tk::real, 3 >& z )
// *****************************************************************************
//! Compute geometry of the face given by three vertices
//! \param[in] x x-coordinates of the three vertices of the triangular face.
//! \param[in] y y-coordinates of the three vertices of the triangular face.
//! \param[in] z z-coordinates of the three vertices of the triangular face.
//! \return Face geometry information. This includes face area, unit normal
//!   pointing outward of the element to the left of the face, and face
//!   centroid coordinates.
// *****************************************************************************
{
  tk::Fields geoiFace( 1, 7 );

  tk::real xp1, yp1, zp1,
           xp2, yp2, zp2,
           xp3, yp3, zp3,
           ax, ay, az,
           bx, by, bz,
           nx, ny, nz,
           sidea, sideb, sidec,
           semip, farea;

  xp1 = x[0];
  xp2 = x[1];
  xp3 = x[2];

  yp1 = y[0];
  yp2 = y[1];
  yp3 = y[2];

  zp1 = z[0];
  zp2 = z[1];
  zp3 = z[2];

  sidea = sqrt( (xp2-xp1)*(xp2-xp1)
              + (yp2-yp1)*(yp2-yp1)
              + (zp2-zp1)*(zp2-zp1) );

  sideb = sqrt( (xp3-xp2)*(xp3-xp2)
              + (yp3-yp2)*(yp3-yp2)
              + (zp3-zp2)*(zp3-zp2) );

  sidec = sqrt( (xp1-xp3)*(xp1-xp3)
              + (yp1-yp3)*(yp1-yp3)
              + (zp1-zp3)*(zp1-zp3) );

  semip = 0.5 * (sidea + sideb + sidec);

  farea = sqrt( semip
              * (semip-sidea)
              * (semip-sideb)
              * (semip-sidec) );

  geoiFace(0,0,0) = farea;

  // get unit normal to face
  ax = xp2 - xp1;
  ay = yp2 - yp1;
  az = zp2 - zp1;

  bx = xp3 - xp1;
  by = yp3 - yp1;
  bz = zp3 - zp1;

  nx =   ay*bz - az*by;
  ny = -(ax*bz - az*bx);
  nz =   ax*by - ay*bx;

  farea = sqrt(nx*nx + ny*ny + nz*nz);

  geoiFace(0,1,0) = nx/farea;
  geoiFace(0,2,0) = ny/farea;
  geoiFace(0,3,0) = nz/farea;

  // get centroid
  geoiFace(0,4,0) = (xp1+xp2+xp3)/3.0;
  geoiFace(0,5,0) = (yp1+yp2+yp3)/3.0;
  geoiFace(0,6,0) = (zp1+zp2+zp3)/3.0;

  return geoiFace;
}
        
tk::Fields
genGeoElemTet( const std::vector< std::size_t >& inpoel,
               const tk::UnsMesh::Coords& coord )
// *****************************************************************************
//  Generate derived data, which stores the geometry details of tetrahedral
//   elements.
//! \param[in] inpoel Element-node connectivity.
//! \param[in] coord Co-ordinates of nodes in this mesh-chunk.
//! \return Element geometry information. This includes element volume and
//!   element centroid coordinates. Use the following examples to access this
//!   information for element-e.
//!   volume: geoElem(e,0,0),
//!   centroid x-coordinate: geoElem(f,1,0),
//!            y-coordinate: geoElem(f,2,0),
//!            z-coordinate: geoElem(f,3,0).
// *****************************************************************************
{
  // set tetrahedron geometry
  std::size_t nnpe(4);

  Assert( inpoel.size()%nnpe == 0,
          "Size of inpoel must be divisible by nnpe" );

  auto nelem = inpoel.size()/nnpe;

  tk::Fields geoElem( nelem, 4 );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for(std::size_t e=0; e<nelem; ++e)
  {
    // get volume
    const auto A = inpoel[nnpe*e+0];
    const auto B = inpoel[nnpe*e+1];
    const auto C = inpoel[nnpe*e+2];
    const auto D = inpoel[nnpe*e+3];
    std::array< tk::real, 3 > ba{{ x[B]-x[A], y[B]-y[A], z[B]-z[A] }},
                              ca{{ x[C]-x[A], y[C]-y[A], z[C]-z[A] }},
                              da{{ x[D]-x[A], y[D]-y[A], z[D]-z[A] }};

    const auto vole = tk::triple( ba, ca, da ) / 6.0;

    Assert( vole > 0, "Element Jacobian non-positive" );

    geoElem(e,0,0) = vole;

    // get centroid
    geoElem(e,1,0) = (x[A]+x[B]+x[C]+x[D])/4.0;
    geoElem(e,2,0) = (y[A]+y[B]+y[C]+y[D])/4.0;
    geoElem(e,3,0) = (z[A]+z[B]+z[C]+z[D])/4.0;
  }

  return geoElem;
}

} // tk::
