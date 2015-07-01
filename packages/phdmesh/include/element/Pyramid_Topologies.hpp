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

#ifndef phdmesh_Pyramid_Topologies_hpp
#define phdmesh_Pyramid_Topologies_hpp

#include <element/CellTopology.hpp>
#include <element/Triangle_Topologies.hpp>
#include <element/Quadrilateral_Topologies.hpp>

namespace phdmesh {

/** \class Pyramid
 *  \brief Topological traits for 5, 13, 14 node pyramids.
 */
template< unsigned = 5 > struct Pyramid ;

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,   5 > ,
                IndexList< 1 , 2 ,   6 > ,
                IndexList< 2 , 0 ,   7 > ,
                IndexList< 0 , 3 ,   8 > ,
                IndexList< 0 , 4 ,   9 > ,
                IndexList< 1 , 4 ,  10 > ,
                IndexList< 2 , 4 ,  11 > ,
                IndexList< 3 , 4 ,  12 > >::type PyramidEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 ,       5 , 10 ,  9 > ,
                IndexList< 1 , 2 , 4 ,       6 , 11 , 10 > ,
                IndexList< 2 , 3 , 4 ,       7 , 12 , 11 > ,
                IndexList< 3 , 0 , 4 ,       8 ,  9 , 12 > ,
                IndexList< 0 , 3 , 2 , 1 ,   8 ,  7 ,  6 ,  5 ,  13 >
  >::type PyramidFaceNodeMap ;

//----------------------------------------------------------------------

template<> struct Pyramid<5> : public
  CellTopologyTraits< 3 , 5 , 5 ,
                      MakeTypeList< Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  >::type ,
                      PyramidEdgeNodeMap ,
                      MakeTypeList< Triangle<>  ,
                                    Triangle<>  ,
                                    Triangle<>  ,
                                    Triangle<>  ,
                                    Quadrilateral<>  >::type ,
                      PyramidFaceNodeMap >
{ typedef Pyramid<5> base ; };

template<> struct Pyramid<13> : public
  CellTopologyTraits< 3 , 5 , 13 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      PyramidEdgeNodeMap ,
                      MakeTypeList< Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Quadrilateral<8>  >::type ,
                      PyramidFaceNodeMap >
{ typedef Pyramid<5> base ; };

template<> struct Pyramid<14> : public
  CellTopologyTraits< 3 , 5 , 14 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      PyramidEdgeNodeMap ,
                      MakeTypeList< Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Quadrilateral<9>  >::type ,
                      PyramidFaceNodeMap >
{ typedef Pyramid<5> base ; };

template<> const CellTopology * cell_topology< Pyramid<> >();
template<> const CellTopology * cell_topology< Pyramid<13> >();
template<> const CellTopology * cell_topology< Pyramid<14> >();

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

