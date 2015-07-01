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

#ifndef phdmesh_Hexahedron_Topologies_hpp
#define phdmesh_Hexahedron_Topologies_hpp

#include <element/CellTopology.hpp>
#include <element/Quadrilateral_Topologies.hpp>

namespace phdmesh {

/** \class Hexahedron
 *  \brief Topological traits for 8, 20, and 27 node hexahedrons.
 */
template< unsigned = 8 > struct Hexahedron ;

//----------------------------------------------------------------------
/*--------------------------------------------------------------------*/
/**
 *  Linear 8-Node Hexahedron Nodes
 *
 *         7                    6
 *          o------------------o
 *         /|                 /|
 *        / |                / |
 *       /  |               /  |
 *      /   |              /   |
 *     /    |             /    |
 *    /     |            /     |
 * 4 /      |         5 /      |
 *  o------------------o       |
 *  |       |          |       |
 *  |     3 o----------|-------o 2
 *  |      /           |      /
 *  |     /            |     /
 *  |    /             |    /
 *  |   /              |   /
 *  |  /               |  /
 *  | /                | /
 *  |/                 |/
 *  o------------------o
 * 0                    1
 *
 *--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/**
 *  Quadratic 20-Node Hexahedron Nodes
 *
 *          7         18         6
 *           o--------o---------o
 *          /|                 /|
 *         / |                / |
 *        /  |               /  |
 *     19o   |            17o   |
 *      /  15o             /    o14
 *     /     |            /     |
 *  4 /      | 16        /      |
 *   o---------o--------o 5     |
 *   |       |       10 |       |
 *   |     3 o-------o--|-------o 2
 *   |      /           |      /
 *   |     /            |     /
 * 12o    /             o13  /
 *   |   o11            |   o9
 *   |  /               |  /
 *   | /                | /
 *   |/                 |/
 *   o---------o--------o
 *  0          8         1
 *
 *--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/**
 *  Quadratic 27-Node Hexahedron Nodes
 *
 *           x--------x---------x
 *          /|                 /|
 *         / |                / |
 *        /  |   22          /  |
 *       x   |    o         x   |
 *      /    x       o26   /    x     Node #20 is at centroid of element
 *     /     |            /     |
 *    /      |           /      |     "2D surface" containing nodes
 *   x---------x--------x       |      0,1,5,4 has node 25 at center....
 *   | 23o   |          |   o24 |
 *   |       x-------x--|-------x
 *   |      /           |      /
 *   |     /  25        |     /
 *   x    /    o        x    /
 *   |   x        o21   |   x
 *   |  /               |  /
 *   | /                | /
 *   |/                 |/
 *   x---------x--------x
 *
 *--------------------------------------------------------------------*/
//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,   8 > ,
                IndexList< 1 , 2 ,   9 > ,
                IndexList< 2 , 3 ,  10 > ,
                IndexList< 3 , 0 ,  11 > ,
                IndexList< 4 , 5 ,  16 > ,
                IndexList< 5 , 6 ,  17 > ,
                IndexList< 6 , 7 ,  18 > ,
                IndexList< 7 , 4 ,  19 > ,
                IndexList< 0 , 4 ,  12 > ,
                IndexList< 1 , 5 ,  13 > ,
                IndexList< 2 , 6 ,  14 > ,
                IndexList< 3 , 7 ,  15 > >::type
  HexahedronEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0, 1, 5, 4,   8, 13, 16, 12,   25 > ,
                IndexList< 1, 2, 6, 5,   9, 14, 17, 13,   24 > ,
                IndexList< 2, 3, 7, 6,  10, 15, 18, 14,   26 > ,
                IndexList< 0, 4, 7, 3,  12, 19, 15, 11,   23 > ,
                IndexList< 0, 3, 2, 1,  11, 10,  9,  8,   21 > ,
                IndexList< 4, 5, 6, 7,  16, 17, 18, 19,   22 > >::type
  HexahedronFaceNodeMap ;

//----------------------------------------------------------------------

template<> struct Hexahedron<8> : public
  CellTopologyTraits< 3 , 8 , 8 ,
                      MakeTypeList< Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  >::type ,
                      HexahedronEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<>  , 
                                    Quadrilateral<>  ,
                                    Quadrilateral<>  ,
                                    Quadrilateral<>  ,
                                    Quadrilateral<>  ,
                                    Quadrilateral<>  >::type ,
                      HexahedronFaceNodeMap >
{
  typedef Hexahedron<8> base ;
};

template<> struct Hexahedron<20> : public
  CellTopologyTraits< 3 , 8 , 20 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      HexahedronEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<8>  , 
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  >::type ,
                      HexahedronFaceNodeMap >
{
  typedef Hexahedron<8> base ;
};

template<> struct Hexahedron<27> : public
  CellTopologyTraits< 3 , 8 , 27 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      HexahedronEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<9>  , 
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  >::type ,
                      HexahedronFaceNodeMap >
{
  typedef Hexahedron<8> base ;
};

template<> const CellTopology * cell_topology< Hexahedron<8> >();
template<> const CellTopology * cell_topology< Hexahedron<20> >();
template<> const CellTopology * cell_topology< Hexahedron<27> >();

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

