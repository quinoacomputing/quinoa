/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_BOXMESHFIXTURE_HPP
#define KOKKOS_BOXMESHFIXTURE_HPP

#include <stdexcept>
#include <KokkosArray_MDArray.hpp>
#include <KokkosArray_MultiVector.hpp>

//  construct a structured, rectangular prism mesh of Hex elements,
//  with dimensions given by elems_x, elems_y, elems_z

template < class IndexArray , class ScalarArray >
class MeshFixture {
public:
  IndexArray   elem_node_ids ;
  IndexArray   node_elem_offset ;
  IndexArray   node_elem_ids ;
  ScalarArray  node_coords ;
};

template < typename Scalar , class Device >
class BoxMeshFixture {
public:
  enum { ELEMENT_NODE_COUNT = 8 };

  typedef typename Device::size_type  index_type ;
  typedef Scalar                      scalar_type ;

  typedef KokkosArray::MDArray< index_type ,  Device > index_array_d ;
  typedef KokkosArray::MDArray< scalar_type , Device > scalar_array_d ;

  typedef typename index_array_d  ::HostMirror  index_array_h ;
  typedef typename scalar_array_d ::HostMirror  scalar_array_h ;

  MeshFixture< index_array_d , scalar_array_d > d_mesh ;
  MeshFixture< index_array_h , scalar_array_h > h_mesh ;

  const index_type elem_count_x ;
  const index_type elem_count_y ;
  const index_type elem_count_z ;
  const index_type elem_count ;
  const index_type node_count_x ;
  const index_type node_count_y ;
  const index_type node_count_z ;
  const index_type node_count ;

  //        7 -------------- 6
  //       /|               /|
  //      / |              / |
  //     /  |             /  |
  //    /   |            /   |
  //   4 -------------- 5    |
  //   |    |           |    |
  //   |    |           |    |
  //   |    3 --------- | -- 2
  //   |   /            |   /
  //   |  /             |  /
  //   | /              | /
  //   |/               |/
  //   0 -------------- 1
  //
  //   Coordinate system:
  //     X = grid-i
  //     Y = grid-j
  //     Z = grid-k
  //
  //   Z   Y
  //   |  /
  //   | /
  //   |/
  //   O-----X
  //

  index_type node_id( const index_type ix , const index_type iy , const index_type iz ) const
  { return ix + node_count_x * ( iy + node_count_y * iz ); }

  index_type elem_id( const index_type ix , const index_type iy , const index_type iz ) const
  { return ix + elem_count_x * ( iy + elem_count_y * iz ); }

  void node_grid( index_type id , index_type & ix , index_type & iy , index_type & iz ) const
  {
    ix = id % node_count_x ; id /= node_count_x ;
    iy = id % node_count_y ; id /= node_count_y ;
    iz = id ;
  }

  void elem_grid( index_type id , index_type & ix , index_type & iy , index_type & iz ) const
  {
    ix = id % elem_count_x ; id /= elem_count_x ;
    iy = id % elem_count_y ; id /= elem_count_y ;
    iz = id ;
  }

  void verify_connectivity_and_coordinates() const
  {
    for ( index_type node_index = 0 ; node_index < node_count; ++node_index ) {
      for ( index_type j = h_mesh.node_elem_offset( node_index ) ;
                j < h_mesh.node_elem_offset( node_index + 1 ) ; ++j ) {
        const index_type elem_index = h_mesh.node_elem_ids(j,0);
        const index_type node_local = h_mesh.node_elem_ids(j,1);
        const index_type en_id = h_mesh.elem_node_ids( elem_index , node_local );
        if ( node_index != en_id ) {
          throw std::runtime_error( std::string("node_elem_ids mapping failure") );
        }
      }
    }

    const index_type elem_node_local_coord[ ELEMENT_NODE_COUNT ][3] =
      { { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 1 , 1 , 0 } , { 0 , 1 , 0 } ,
        { 0 , 0 , 1 } , { 1 , 0 , 1 } , { 1 , 1 , 1 } , { 0 , 1 , 1 } };

    for ( index_type elem_index = 0 ; elem_index < elem_count; ++elem_index ) {

      index_type g[3] ;

      elem_grid( elem_index , g[0] , g[1] , g[2] );

      for ( index_type nn = 0 ; nn < ELEMENT_NODE_COUNT ; ++nn ) {
        const index_type node_index = h_mesh.elem_node_ids( elem_index , nn );

        for ( index_type nc = 0 ; nc < 3 ; ++nc ) {
          const index_type en_coord = g[nc] + elem_node_local_coord[ nn ][ nc ];
          const index_type n_coord  = (int)   h_mesh.node_coords( node_index , nc );

          if ( en_coord != n_coord ) {
            throw std::runtime_error( std::string("elem_node_coord mapping failure") );
          }
        }
      }
    }
  }

private:

  void populate_node_element()
  {
    index_array_h node_elem_count = KokkosArray::create_mdarray< index_array_h >( node_count );

    for(index_type i = 0; i < node_count ; i++) {
      node_elem_count( i ) = 0 ;
    }

    for ( index_type i = 0 ; i < elem_count ; ++i ) {
      for ( index_type n = 0 ; n < ELEMENT_NODE_COUNT  ; ++n ) {
        ++node_elem_count( h_mesh.elem_node_ids(i,n) );
      }
    }

    h_mesh.node_elem_offset(0) = 0 ;
    for(index_type i = 0; i < node_count ; ++i ){
      h_mesh.node_elem_offset( i + 1 ) = h_mesh.node_elem_offset(i) + node_elem_count(i);
      node_elem_count( i ) = 0 ;
    }

    // Looping in element order insures the list of elements
    // is sorted by element index.

    for ( index_type i = 0 ; i < elem_count ; ++i ) {
      for ( index_type n = 0 ; n < ELEMENT_NODE_COUNT ; ++n ) {
        const index_type nid = h_mesh.elem_node_ids(i, n);
        const index_type j = h_mesh.node_elem_offset(nid) + node_elem_count(nid);

        h_mesh.node_elem_ids( j , 0 ) = i ;
        h_mesh.node_elem_ids( j , 1 ) = n ;

        ++node_elem_count( nid );
      }
    }
  }

public:

  BoxMeshFixture( const index_type elems_x ,
                  const index_type elems_y ,
                  const index_type elems_z )
  : d_mesh() , h_mesh()
  , elem_count_x( elems_x )
  , elem_count_y( elems_y )
  , elem_count_z( elems_z )
  , elem_count(   elems_z * elems_y * elems_x )
  , node_count_x( elems_x + 1 )
  , node_count_y( elems_y + 1 )
  , node_count_z( elems_z + 1 )
  , node_count( ( elems_z + 1 ) * ( elems_y + 1 ) * ( elems_x + 1 ) )
  {
    const index_type elem_node_local_coord[ ELEMENT_NODE_COUNT ][3] =
      { { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 1 , 1 , 0 } , { 0 , 1 , 0 } ,
        { 0 , 0 , 1 } , { 1 , 0 , 1 } , { 1 , 1 , 1 } , { 0 , 1 , 1 } };

    const index_type count_node_elem = elem_count * ELEMENT_NODE_COUNT ;

    d_mesh.node_coords      = KokkosArray::create_mdarray< scalar_array_d >(node_count, 3);
    d_mesh.elem_node_ids    = KokkosArray::create_mdarray< index_array_d >( elem_count,  ELEMENT_NODE_COUNT );
    d_mesh.node_elem_offset = KokkosArray::create_mdarray< index_array_d >( node_count + 1 );
    d_mesh.node_elem_ids    = KokkosArray::create_mdarray< index_array_d >( count_node_elem , 2 );

    h_mesh.node_coords      = KokkosArray::create_mirror( d_mesh.node_coords );
    h_mesh.elem_node_ids    = KokkosArray::create_mirror( d_mesh.elem_node_ids );
    h_mesh.node_elem_offset = KokkosArray::create_mirror( d_mesh.node_elem_offset );
    h_mesh.node_elem_ids    = KokkosArray::create_mirror( d_mesh.node_elem_ids );

    // Initialize node coordinates of grid.

    for ( index_type node_index = 0 ; node_index < node_count ; ++node_index ) {
      index_type ig , jg , kg ;
      node_grid( node_index , ig , jg , kg );
      h_mesh.node_coords(node_index,0) = ig ;
      h_mesh.node_coords(node_index,1) = jg ;
      h_mesh.node_coords(node_index,2) = kg ;
    }

    // Initialize element-node connectivity:

    for ( index_type elem_index = 0 ; elem_index < elem_count ; ++elem_index ) {
      index_type ig , jg , kg ;
      elem_grid( elem_index , ig , jg , kg );

      for ( index_type nn = 0 ; nn < ELEMENT_NODE_COUNT ; ++nn ) {
        h_mesh.elem_node_ids( elem_index , nn ) =
          node_id( ig + elem_node_local_coord[nn][0] ,
                   jg + elem_node_local_coord[nn][1] ,
                   kg + elem_node_local_coord[nn][2] );
      }
    }

    populate_node_element();

    verify_connectivity_and_coordinates();

    KokkosArray::deep_copy( d_mesh.node_coords ,      h_mesh.node_coords );
    KokkosArray::deep_copy( d_mesh.elem_node_ids ,    h_mesh.elem_node_ids );
    KokkosArray::deep_copy( d_mesh.node_elem_offset , h_mesh.node_elem_offset );
    KokkosArray::deep_copy( d_mesh.node_elem_ids ,    h_mesh.node_elem_ids );
  }


  template< typename ValueType >
  void init_dirichlet_z( KokkosArray::MultiVector< index_type, Device > & node_flag ,
                         KokkosArray::MultiVector< ValueType , Device > & node_value ) const
  {
    typedef KokkosArray::MultiVector< index_type, Device > index_vector_d ;
    typedef KokkosArray::MultiVector< ValueType , Device > value_vector_d ;

    typedef typename value_vector_d::HostMirror  value_vector_h ;
    typedef typename index_vector_d::HostMirror  index_vector_h ;

    node_flag  = KokkosArray::create_multivector< index_vector_d >( node_count );
    node_value = KokkosArray::create_multivector< value_vector_d >( node_count );

    index_vector_h flag_h  = KokkosArray::create_mirror( node_flag );
    value_vector_h value_h = KokkosArray::create_mirror( node_value );

    index_type ig , jg , kg ;
    for ( index_type i = 0 ; i < node_count ; ++i ) {
      node_grid( i , ig , jg , kg );
      if ( kg == 0 || ( kg + 1 ) == node_count_z ) {
        flag_h(i)  = 1 ;
        value_h(i) = kg ;
      }
    }

    KokkosArray::deep_copy( node_flag ,  flag_h );
    KokkosArray::deep_copy( node_value , value_h );
  }
};

#endif /* #ifndef KOKKOS_BOXMESHFIXTURE_HPP */

