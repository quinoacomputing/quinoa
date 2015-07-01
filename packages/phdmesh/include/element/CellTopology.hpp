/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#ifndef phdmesh_CellTopology_hpp
#define phdmesh_CellTopology_hpp

#include <iosfwd>
#include <util/TypeList.hpp>
#include <util/IndexList.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** \struct CellTopology
 *  \brief Simple 'C' struct of cell topology attributes.
 *
 *  The topology may be extended such that the number of nodes
 *  (subcells of dimension zero) is greater than the number of
 *  vertices.  In this case the vertices must be ordered first.
 *
 *  Given a CellTopology object 'top' the attributes are accessed as follows.
 *
 *  - top.base                       ; Base (non-extended) topology
 *  - top.name
 *  - top.key
 *  - top.dimension
 *  - top.vertex_count
 *  - top.node_count
 *  - top.edge_count
 *  - top.side_count
 *
 *  - top.subcell_homogeneity        ; subcell[Dim] are homogeneous
 *  - top.subcell_count[Dim]         ; number of subcells of dimension Dim
 *  - top.subcell[Dim][Ord].topology ; topology of subcell
 *  - top.subcell[Dim][Ord].node[I]  ; node number of subcell's node I
 *
 *  For 2D or 3D cells:
 *  - top.side[Ord].topology         ; subcell[dimension-1][Ord].topology
 *  - top.side[Ord].node[I]          ; subcell[dimension-1][Ord].node[I]
 *  - top.edge[Ord].topology         ; subcell[1][Ord].topology
 *  - top.edge[Ord].node[I]          ; subcell[1][Ord].node[I]
 *
 * Nodes, edges, and sides are subcells with a particular dimension.
 *  - node has Dim == 0
 *  - edge has Dim == 1 when top.dimension >= 2
 *  - side has Dim == top.dimension - 1.
 */
struct CellTopology {
  /** \brief Base, a.k.a. not-extended, version of this topology
   *         where vertex_count == edge_count.
   */
  const CellTopology * base ;

  /** \brief Intuitive name for this topology */
  const char * name ;

  /** \brief Unique key for this topology */
  unsigned key ;

  /** \brief Topological dimension */
  unsigned dimension ;

  /** \brief Number of vertices.
   *  For non-extended topology the number of Cell^0 subcells.
   */
  unsigned vertex_count ;
  unsigned node_count ;
  unsigned edge_count ;
  unsigned side_count ;

  /** \brief Flag if the subcells of a given dimension are homogeneous */
  unsigned subcell_homogeneity ;

  /** \brief Number of subcells of each dimension. */
  unsigned subcell_count[4] ;

  /** \brief Subcell information.
   *  - required: 0 <= Ord <= subcell_count[Dim]
   *  - required: 0 <= J   <  subcell[Dim][Ord]->subcell_count[0]
   *  - subcell[Dim][Ord].topology
   *  - subcell[Dim][Ord].node[J]
   */
  struct Subcell {
    /** \brief Subcell topology */
    const CellTopology * topology ;

    /** \brief Subcell indexing of Cell^0 with respect to parent cell. */
    const unsigned * node ;
  };

  /** \brief Subcells of each dimension */
  const struct Subcell * subcell[4] ;

  const struct Subcell * side ;
  const struct Subcell * edge ;
};

template< class Traits >
const CellTopology * cell_topology();

std::ostream & operator << ( std::ostream & , const CellTopology & );

//----------------------------------------------------------------------
/** \struct CellTopology Traits
 *  \brief  Cell topology attributes as compile-time traits.
 *
 * Given a CellTopologyTraits 'Top' the attributes are accessed as follows.
 *  - Top::key
 *  - Top::dimension
 *  - Top::vertex_count
 *  - top::node_count
 *  - top::edge_count
 *  - top::side_count
 *  - top::subcell_homogeneity
 *
 *  - Top::subcell<Dim>::count        ;  number of subcells of dimension Dim
 *  - Top::subcell<Dim,Ord>::topology ;  topology traits of subcell
 *  - Top::subcell<Dim,Ord,I>::node   ;  node number of subcell's node I
 *
 *  - Top::side<Ord>::topology        ; subcell<dimension-1,Ord>::topology
 *  - Top::side<Ord,I>::node          ; subcell<dimension-1,Ord,I>::.node
 *  - Top::edge<Ord>::topology        ; subcell<1,Ord>::topology
 *  - Top::edge<Ord,I>::node          ; subcell<1,Ord,I>::.node
 */
template< unsigned Dimension ,
          unsigned Number_Vertex ,
          unsigned Number_Node ,
          class    EdgeList = TypeListEnd ,
          class    EdgeMaps = TypeListEnd ,
          class    FaceList = TypeListEnd ,
          class    FaceMaps = TypeListEnd >
struct CellTopologyTraits ;

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Implementation details to follow.

namespace phdmesh {

template< class CellTop , class CellMap , unsigned Index , bool Good >
struct SubcellNodeIndex ;

template< class CellTop , class CellMap , unsigned Index >
struct SubcellNodeIndex< CellTop , CellMap , Index , false >
{ enum { value = -1 }; };

template< class CellTop , class CellMap , unsigned Index >
struct SubcellNodeIndex< CellTop , CellMap , Index , true >
{
  enum { value = Index < CellTop::template subcell<0>::count
               ? IndexListAt< CellMap , Index >::value : -1 };
};

//----------------------------------------------------------------------

template< unsigned SubcellDim , unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned Dimension ,
          unsigned Number_Vertex , unsigned Number_Node ,
          class EdgeList , class EdgeMaps ,
          class FaceList , class FaceMaps >
struct SubcellTopologyTraits ;

//----------------------------------------------------------------------
// Self-subcell reference

template< unsigned NodeIndex ,
          unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps >
struct SubcellTopologyTraits<0,0,NodeIndex, 0,NV,NN,EList,EMaps,FList,FMaps>
{
  typedef CellTopologyTraits<0,NV,NN,EList,EMaps,FList,FMaps> topology ;
  enum { count = 1 };
  enum { node = NodeIndex < NN ? (int) NodeIndex : -1 };
};

template< unsigned NodeIndex ,
          unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps >
struct SubcellTopologyTraits<1,0,NodeIndex, 1,NV,NN,EList,EMaps,FList,FMaps>
{
  typedef CellTopologyTraits<1,NV,NN,EList,EMaps,FList,FMaps> topology ;
  enum { count = 1 };
  enum { node = NodeIndex < NN ? (int) NodeIndex : -1 };
};

template< unsigned NodeIndex ,
          unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps >
struct SubcellTopologyTraits<2,0,NodeIndex, 2,NV,NN,EList,EMaps,FList,FMaps>
{
  typedef CellTopologyTraits<2,NV,NN,EList,EMaps,FList,FMaps> topology ;
  enum { count = 1 };
  enum { node = NodeIndex < NN ? (int) NodeIndex : -1 };
};

template< unsigned NodeIndex ,
          unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps >
struct SubcellTopologyTraits<3,0,NodeIndex, 3,NV,NN,EList,EMaps,FList,FMaps>
{
  typedef CellTopologyTraits<3,NV,NN,EList,EMaps,FList,FMaps> topology ;
  enum { count = 1 };
  enum { node = NodeIndex < NN ? (int) NodeIndex : -1 };
};

//----------------------------------------------------------------------
// Node-subcell reference:

template< unsigned SubcellOrd ,
          unsigned D , unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps >
struct SubcellTopologyTraits<0,SubcellOrd,0, D,NV,NN,EList,EMaps,FList,FMaps>
{
  typedef CellTopologyTraits<0,0,0> topology ;
  enum { count = NN };
  enum { node = SubcellOrd < NN ? (int) SubcellOrd : -1 };
};

// Edge-subcell reference:

template< unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned D , unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps >
struct SubcellTopologyTraits<1,SubcellOrd,NodeIndex,
                             D,NV,NN,EList,EMaps,FList,FMaps>
{
private:
  typedef typename TypeListAt<EMaps,SubcellOrd>::type node_map ;
public:

  typedef typename TypeListAt<EList,SubcellOrd>::type topology ;

  enum { count = TypeListLength<EList>::value };

  enum { node = SubcellNodeIndex< topology , node_map , NodeIndex ,
                                  SubcellOrd < count >::value };
};

// Face-subcell reference:

template< unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned D , unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps >
struct SubcellTopologyTraits< 2, SubcellOrd, NodeIndex,
                              D,NV,NN,EList,EMaps,FList,FMaps>
{
private:
  typedef typename TypeListAt<FMaps,SubcellOrd>::type node_map ;
public:

  typedef typename TypeListAt<FList,SubcellOrd>::type topology ;

  enum { count = TypeListLength<FList>::value };

  enum { node = SubcellNodeIndex< topology , node_map , NodeIndex ,
                                  SubcellOrd < count >::value };
};

//----------------------------------------------------------------------
// Only partially specialized subcell references are valid.

template< unsigned SubcellDim , unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned Dimension ,
          unsigned Number_Vertex , unsigned Number_Node ,
          class EdgeList , class EdgeMaps ,
          class FaceList , class FaceMaps >
struct SubcellTopologyTraits
{
  typedef void topology ;
  enum { count = 0 };
  enum { node = -1 };
};

//----------------------------------------------------------------------

template< class ListType > struct TypeListHomogeneous ;

template<>
struct TypeListHomogeneous<TypeListEnd> {
  enum { value = true };
};

template< class T >
struct TypeListHomogeneous< TypeList<T,TypeListEnd> > {
  enum { value = true };
};

template< class T , class Tail >
struct TypeListHomogeneous< TypeList< T, TypeList< T , Tail > > > {
  enum { value = TypeListHomogeneous< TypeList<T,Tail> >::value };
};

template< class ListType >
struct TypeListHomogeneous 
{
  enum { value = false };
};

//----------------------------------------------------------------------

template< unsigned Dimension , unsigned Number_Vertex , unsigned Number_Node ,
          class EdgeList , class EdgeMaps ,
          class FaceList , class FaceMaps >
struct CellTopologyTraits {

  typedef CellTopologyTraits< Dimension, Number_Vertex, Number_Node,
                              EdgeList, EdgeMaps, FaceList, FaceMaps > Traits ;

  enum { dimension    = Dimension ,
         vertex_count = Number_Vertex ,
         node_count   = Number_Node ,
         edge_count   = TypeListLength<EdgeList>::value ,
         face_count   = TypeListLength<FaceList>::value ,
         side_count   = Dimension == 3 ? face_count : (
                        Dimension == 2 ? edge_count : 0 ),
         key          = ( dimension    << 28 ) |
                        ( face_count   << 24 ) |
                        ( edge_count   << 20 ) |
                        ( vertex_count << 16 ) |
                        ( node_count ) };

  enum { subcell_homogeneity = TypeListHomogeneous<EdgeList>::value &&
                               TypeListHomogeneous<FaceList>::value };

  template< unsigned Dim, unsigned Ord = 0, unsigned J = 0 >
  struct subcell :
    public SubcellTopologyTraits< Dim , Ord , J ,
                                  dimension , vertex_count , node_count ,
                                  EdgeList , EdgeMaps ,
                                  FaceList , FaceMaps > {};

  template< unsigned Ord = 0 , unsigned J = 0 >
  struct side :
    public SubcellTopologyTraits< ( 1 < dimension ? dimension - 1 : 4 ) ,
                                  Ord , J ,
                                  dimension , vertex_count , node_count ,
                                  EdgeList , EdgeMaps ,
                                  FaceList , FaceMaps > {};

  template< unsigned Ord = 0 , unsigned J = 0 >
  struct edge :
    public SubcellTopologyTraits< ( 1 < dimension ? 1 : 4 ) , Ord , J ,
                                  dimension , vertex_count , node_count ,
                                  EdgeList , EdgeMaps ,
                                  TypeListEnd , TypeListEnd > {};

private:

  enum { nedge_map = TypeListLength<EdgeMaps>::value ,
         nface_map = TypeListLength<FaceMaps>::value };

  enum { OK_edge = StaticAssert< (int) edge_count == (int) nedge_map >::OK };
  enum { OK_face = StaticAssert< (int) face_count == (int) nface_map >::OK };
  enum { OK_key  = StaticAssert< ( 0 == Dimension     >> 4 ) &&
                                 ( 0 == face_count   >> 4 ) &&
                                 ( 0 == edge_count   >> 4 ) &&
                                 ( 0 == Number_Vertex >> 4 ) &&
                                 ( 0 == Number_Node   >> 16 ) >::OK };
};

} // namespace phdmesh

#endif

