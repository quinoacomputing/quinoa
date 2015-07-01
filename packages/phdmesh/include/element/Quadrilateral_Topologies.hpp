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

#ifndef phdmesh_Quadrilateral_Topologies_hpp
#define phdmesh_Quadrilateral_Topologies_hpp

#include <element/CellTopology.hpp>
#include <element/Basic_Topologies.hpp>

namespace phdmesh {

template< unsigned = 4 > struct Quadrilateral ;
template< unsigned = 4 > struct ShellQuadrilateral ;

//----------------------------------------------------------------------
// Conventional numbering quadrilateral with up to nine-nodes
//
//                 Edge #2
//
//            3        6        2
//             o-------o-------o
//             |               |
//             |               |
//             |       8       |
//  Edge #3  7 o       o       o 5  Edge #1
//             |               |
//             |               |
//             |               |
//             o-------o-------o
//            0        4        1
//
//                  Edge #0
//
// Conformal node numbering with up to 25 nodes:
//
//    3    14    6   13     2
//     o----*----o----*----o
//     |         |         |
//     |   24    |    23   |
//   15*    *    *19  *    *12
//     |         |         |
//     |        8|    18   |
//   7 o----*----o----*----o 5
//     |   20    |         |
//     |         |         |
//   16*    *  17*    *    *11
//     |   21    |   22    |
//     |         |         |
//     o----*----o----*----o
//    0     9    4   10     1
//
//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,  4 ,   9 , 10 > ,
                IndexList< 1 , 2 ,  5 ,  11 , 12 > ,
                IndexList< 2 , 3 ,  6 ,  13 , 14 > ,
                IndexList< 3 , 0 ,  7 ,  15 , 16 > >::type
  QuadrilateralEdgeNodeMap ;

//----------------------------------------------------------------------

template<> struct Quadrilateral<4> : public
  CellTopologyTraits< 2 , 4 , 4 ,
                            MakeTypeList< Line<>  ,
                                          Line<>  ,
                                          Line<>  ,
                                          Line<>  >::type ,
                            QuadrilateralEdgeNodeMap >
{ typedef Quadrilateral<4> base ; };

template<> struct Quadrilateral<8> : public
  CellTopologyTraits< 2 , 4 , 8 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap >
{ typedef Quadrilateral<4> base ; };

template<> struct Quadrilateral<9> : public
  CellTopologyTraits< 2 , 4 , 9 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap >
{ typedef Quadrilateral<4> base ; };

template<> const CellTopology * cell_topology< Quadrilateral<> >();
template<> const CellTopology * cell_topology< Quadrilateral<8> >();
template<> const CellTopology * cell_topology< Quadrilateral<9> >();

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0,1,2,3,  4,5,6,7,  8 > ,
                IndexList< 0,3,2,1,  7,6,5,4,  8 > >::type
    ShellQuadrilateralFaceNodeMap ;

template<> struct ShellQuadrilateral<4> : public
  CellTopologyTraits< 3 , 4 , 4 ,
                            MakeTypeList< Line<>  ,
                                          Line<>  ,
                                          Line<>  ,
                                          Line<>  >::type ,
                            QuadrilateralEdgeNodeMap ,
                            MakeTypeList< Quadrilateral<>  ,
                                          Quadrilateral<>  >::type ,
                            ShellQuadrilateralFaceNodeMap >
{ typedef ShellQuadrilateral<4> base ; };

template<> struct ShellQuadrilateral<8> : public
  CellTopologyTraits< 3 , 4 , 8 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap ,
                            MakeTypeList< Quadrilateral<8>  ,
                                          Quadrilateral<8>  >::type ,
                            ShellQuadrilateralFaceNodeMap >
{ typedef ShellQuadrilateral<4> base ; };

template<> struct ShellQuadrilateral<9> : public
  CellTopologyTraits< 3 , 4 , 9 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap ,
                            MakeTypeList< Quadrilateral<9>  ,
                                          Quadrilateral<9>  >::type ,
                            ShellQuadrilateralFaceNodeMap >
{ typedef ShellQuadrilateral<4> base ; };

template<> const CellTopology * cell_topology< ShellQuadrilateral<> >();
template<> const CellTopology * cell_topology< ShellQuadrilateral<8> >();
template<> const CellTopology * cell_topology< ShellQuadrilateral<9> >();

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

