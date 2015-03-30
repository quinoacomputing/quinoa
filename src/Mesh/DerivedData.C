//******************************************************************************
/*!
  \file      src/Mesh/DerivedData.C
  \author    J. Bakosi
  \date      Mon 30 Mar 2015 06:49:03 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Generate data structures derived from unstructured mesh
  \details   Generate data structures derived from the connectivity information
     of an unstructured mesh.
*/
//******************************************************************************

#include <algorithm>
#include <set>
#include <map>

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

  // Return (move out) linked lists
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

  // Return (move out) linked lists
  return std::make_pair( std::move(psup1), std::move(psup2) );
}

std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
genEdsup( const std::vector< int >& inpoel,
          std::size_t nnpe,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup )
//******************************************************************************
//  Generate derived data structure, edges surrounding points
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< int > inpoel { 12, 14,  9, 11,
//!                                 10, 14, 13, 12 };
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
//!          use edge with point ids p < edsup.first[i]
//!   \endcode
//!     To find out the number of points, _npoin_, the mesh connectivity,
//!     _inpoel_, can be queried:
//!   \code{.cpp}
//!     auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
//!     Assert( *minmax.first == 0, "node ids should start from zero" );
//!     auto npoin = static_cast< std::size_t >( *minmax.second + 1 );
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see tk::genInpoed for similar data that sometimes may be more advantageous
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genEdsup() on empty container" );
  Assert( nnpe > 0, "Attempt to call genEdsup() with zero nodes per element" );
  Assert( !esup.first.empty(), "Attempt to call genEdsup() with empty esup1" );
  Assert( !esup.second.empty(), "Attempt to call genEdsup() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = static_cast< std::size_t >( *minmax.second + 1 );

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
        auto q = static_cast< std::size_t >( inpoel[ esup1[i] * nnpe + n ] );
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
  edsup2.push_back( edsup2.back() );

  // Return (move out) linked lists
  return std::make_pair( std::move(edsup1), std::move(edsup2) );
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
//!   access to particular edges. Accessing a unique edge i using the edge
//!   connectivity data structure, inpoed, generated here can be accomplished by
//!   \code{.cpp}
//!      edge i point id p = inpoed[i*2];
//!      edge i point id q = inpoed[i*2+1];
//!   \endcode
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see tk::genEdsup for similar data that sometimes may be more advantageous
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genInpoed() on empty container" );
  Assert( nnpe > 0, "Attempt to call genInpoed() with zero nodes per element" );
  Assert( !esup.first.empty(), "Attempt to call genInpoed() with empty esup1" );
  Assert( !esup.second.empty(),
          "Attempt to call genInpoed() with empty esup2" );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = static_cast< std::size_t >( *minmax.second + 1 );

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
        auto q = static_cast< std::size_t >( inpoel[ esup1[i] * nnpe + n ] );
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

  // Return (move out) linked lists
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
  Assert( !esupel.first.empty(),
          "Attempt to call genEsupel() with empty esupel1" );
  Assert( !esupel.second.empty(),
          "Attempt to call genEsupel() with empty esupel2" );

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

  // Return (move out) vector
  return esuel;
}

std::vector< std::size_t >
genInedel( const std::vector< int >& inpoel,
           std::size_t nnpe,
           const std::vector< std::size_t >& inpoed )
//******************************************************************************
//  Generate derived data structure, edges of elements
//! \param[in] inpoel Inteconnectivity of points and elements. These are the
//!   node ids of each element of an unstructured mesh. Example:
//!   \code{.cpp}
//!     std::vector< int > inpoel { 12, 14,  9, 11,
//!                                 10, 14, 13, 12 };
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
//!   edge ids * 2 (as defined by inpoed) of all elements. The edge ids are
//!   multiplied by 2 so that they can be directly used to index the vector
//!   inpoed. Because the derived data structure generated here, inedel, is
//!   intended to be used in conjunction with the linear vector inpoed and not
//!   with the linked lists edsup1 and edsup2, this function takes inpoed as an
//!   argument.
//! \note At first sight, this function seems to work for elements with more
//!   vertices than that of tetrahedra. However, that is not the case since the
//!   algorithm for nnpe > 4 would erronously identify any two combination of
//!   vertices as a valid edge of an element. Since only triangles and
//!   tetrahedra have no internal edges, this algorithm only works for triangle
//!   and tetrahedra element connectivity.
//! \see Lohner, An Introduction to Applied CFD Techniques, Wiley, 2008
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !inpoel.empty(), "Attempt to call genInedel() on empty container" );
  Assert( nnpe > 0, "Attempt to call genInedel() with zero nodes per element" );
  Assert( !inpoed.empty(), "Attempt to call genInedel() with empty inpoed" );

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

  // generate stars from inpoed; every odd is a star center, every even is a
  // spike
  for (std::size_t i=0; i<inpoed.size()/2; ++i)
    star[ inpoed[i*2] ].push_back( inpoed[i*2+1] );

  // store index of star centers in vector; assume non-center points of each
  // star have already been sorted
  std::vector< std::size_t > edsup2( 1, 0 );
  for (auto p : star) edsup2.push_back( edsup2.back() + p.second.size() );
  edsup2.push_back( edsup2.back() );
  star.clear();

  // Second, generate edges of elements

  auto nelem = inpoel.size()/nnpe;

  // map associating elem id with vector of edge ids
  std::map< std::size_t, std::vector< std::size_t > > edges;

  // generate map of elements associated to edge ids
  for (std::size_t e=0; e<nelem; ++e)
    for (std::size_t n=0; n<nnpe; ++n) {
      auto p = static_cast< std::size_t >( inpoel[e*nnpe+n] );
      for (auto i=edsup2[p]+1; i<=edsup2[p+1]; ++i)
         for (std::size_t j=0; j<nnpe; ++j)
            if (inpoed[(i-1)*2+1] == static_cast<std::size_t>(inpoel[e*nnpe+j]))
               edges[e].push_back( (i-1)*2 );
    }

  // linear vector to store the edge ids of all elements
  std::vector< std::size_t > inedel;

  // store edge ids of elements in linear vector
  for (auto e : edges) for (auto p : e.second) inedel.push_back( p );

  // Return (move out) vector
  return inedel;
}

} // tk::
