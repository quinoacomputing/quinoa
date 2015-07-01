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

#ifndef phdmesh_Triangle_Topologies_hpp
#define phdmesh_Triangle_Topologies_hpp

#include <element/CellTopology.hpp>
#include <element/Basic_Topologies.hpp>

namespace phdmesh {

template< unsigned = 3 > struct Triangle ;
template< unsigned = 3 > struct ShellTriangle ;

/*----------------------------------------------------------------------*/
/*  Triangle node numbering, up to 6 nodes:                             */
/*                                                                      */
/*                  2                                                   */
/*                  o                                                   */
/*                 / \                                                  */
/*                /   \                                                 */
/*               /     \                                                */
/* Edge #2    5 o       o 4   Edge #1                                   */
/*             /         \                                              */
/*            /           \                                             */
/*           /             \                                            */
/*          o-------o-------o                                           */
/*         0        3        1                                          */
/*                                                                      */
/*                Edge #0                                               */
/*                                                                      */
/*  Triangle node numbering, up to 15 nodes conformal to 6 node:        */
/*                                                                      */
/*           2                                                          */
/*           o                                                          */
/*          / \                                                         */
/*      10 *   * 9                                                      */
/*        / 14  \                                                       */
/*     5 o---*---o 4                                                    */
/*      / \     / \                                                     */
/*  11 * 12*   *13 * 8                                                  */
/*    /     \ /     \                                                   */
/*   o---*---o---*---o                                                  */
/*  0    6   3   7    1                                                 */
/*                                                                      */
/*  Interior node 12 is on the line between node 0 and node 4           */
/*  Interior node 13 is on the line between node 1 and node 5           */
/*  Interior node 14 is on the line between node 2 and node 3           */
/*                                                                      */
/*----------------------------------------------------------------------*/

typedef
  MakeTypeList< IndexList< 0 , 1 , 3 ,   6 ,  7 > ,
                IndexList< 1 , 2 , 4 ,   8 ,  9 > ,
                IndexList< 2 , 0 , 5 ,  10 , 11 > >::type
    TriangleEdgeNodeMap ;

//------------------------------------------------------------------------

template<> struct Triangle<3> : public
  CellTopologyTraits< 2 , 3 , 3 ,
                            MakeTypeList< Line<>  ,
                                          Line<>  ,
                                          Line<>  >::type ,
                            TriangleEdgeNodeMap >
{ typedef Triangle<3> base ; };

template<> struct Triangle<6> : public
  CellTopologyTraits< 2 , 3 , 6 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            TriangleEdgeNodeMap >
{ typedef Triangle<3> base ; };

template<> const CellTopology * cell_topology< Triangle<> >();
template<> const CellTopology * cell_topology< Triangle<6> >();

//------------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0,1,2,  3,4,5,   6, 7,8,9,10,11,  12,13,14> ,
                IndexList< 0,2,1,  5,4,3,  11,10,9,8, 7, 6,  12,14,13>
    >::type ShellTriangleFaceNodeMap ;

template<> struct ShellTriangle<3> : public
  CellTopologyTraits< 3 , 3 , 3 ,
                            MakeTypeList< Line<>  ,
                                          Line<>  ,
                                          Line<>  >::type ,
                            TriangleEdgeNodeMap ,
                            MakeTypeList< Triangle<>  ,
                                          Triangle<>  >::type ,
                            ShellTriangleFaceNodeMap >
{ typedef ShellTriangle<3> base ; };

template<> struct ShellTriangle<6> : public
  CellTopologyTraits< 3 , 3 , 6 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            TriangleEdgeNodeMap ,
                            MakeTypeList< Triangle<6>  ,
                                          Triangle<6>  >::type ,
                            ShellTriangleFaceNodeMap >
{ typedef ShellTriangle<3> base ; };

template<> const CellTopology * cell_topology< ShellTriangle<> >();
template<> const CellTopology * cell_topology< ShellTriangle<6> >();

} // namespace phdmesh


#endif

