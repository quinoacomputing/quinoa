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

#ifndef phdmesh_Wedge_Topologies_hpp
#define phdmesh_Wedge_Topologies_hpp

#include <element/CellTopology.hpp>
#include <element/Triangle_Topologies.hpp>
#include <element/Quadrilateral_Topologies.hpp>

namespace phdmesh {

template< unsigned = 6 > struct Wedge ;

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,   6 > ,
                IndexList< 1 , 2 ,   7 > ,
                IndexList< 2 , 0 ,   8 > ,
                IndexList< 3 , 4 ,  12 > ,
                IndexList< 4 , 5 ,  13 > ,
                IndexList< 5 , 3 ,  14 > ,
                IndexList< 0 , 3 ,   9 > ,
                IndexList< 1 , 4 ,  10 > ,
                IndexList< 2 , 5 ,  11 >
  >::type WedgeEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 , 3 ,   6 , 10 , 12 ,  9 ,  15 > ,
                IndexList< 1 , 2 , 5 , 4 ,   7 , 11 , 13 , 10 ,  16 > ,
                IndexList< 0 , 3 , 5 , 2 ,   9 , 14 , 11 ,  8 ,  17 > ,
                IndexList< 0 , 2 , 1 ,       8 ,  7 ,  6 > ,
                IndexList< 3 , 4 , 5 ,      12 , 13 , 14 >
  >::type WedgeFaceNodeMap ;

//----------------------------------------------------------------------

template<> struct Wedge<6> : public
  CellTopologyTraits< 3 , 6 , 6 ,
                      MakeTypeList< Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  >::type ,
                      WedgeEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<>  , 
                                    Quadrilateral<>  ,
                                    Quadrilateral<>  ,
                                    Triangle<>  ,
                                    Triangle<>  >::type ,
                      WedgeFaceNodeMap >
{ typedef Wedge<6> base ; };

template<> struct Wedge<15> : public
  CellTopologyTraits< 3 , 6 , 15 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      WedgeEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<8>  , 
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  >::type ,
                      WedgeFaceNodeMap >
{ typedef Wedge<6> base ; };

template<> struct Wedge<18> : public
  CellTopologyTraits< 3 , 6 , 18 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      WedgeEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<9>  , 
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  >::type ,
                      WedgeFaceNodeMap >
{ typedef Wedge<6> base ; };

template<> const CellTopology * cell_topology< Wedge<> >();
template<> const CellTopology * cell_topology< Wedge<15> >();
template<> const CellTopology * cell_topology< Wedge<18> >();

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

