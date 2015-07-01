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

#ifndef phdmesh_Tetrahedron_Topologies_hpp
#define phdmesh_Tetrahedron_Topologies_hpp

#include <element/CellTopology.hpp>
#include <element/Triangle_Topologies.hpp>

namespace phdmesh {

template< unsigned = 4 > struct Tetrahedron ;

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 > ,
                IndexList< 1 , 2 , 5 > ,
                IndexList< 2 , 3 , 6 > ,
                IndexList< 0 , 3 , 7 > ,
                IndexList< 1 , 3 , 8 > ,
                IndexList< 2 , 3 , 9 > >::type TetrahedronEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 3 ,   4 , 8 , 7 > ,
                IndexList< 1 , 2 , 3 ,   5 , 9 , 8 > ,
                IndexList< 0 , 3 , 2 ,   7 , 9 , 6 > ,
                IndexList< 0 , 2 , 1 ,   6 , 5 , 4 >
  >::type TetrahedronSideNodeMap ;

template<> struct Tetrahedron<4> : public
  CellTopologyTraits< 3 , 4 , 4 ,
                      MakeTypeList< Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  ,
                                    Line<>  >::type ,
                      TetrahedronEdgeNodeMap ,
                      MakeTypeList< Triangle<>  ,
                                    Triangle<>  ,
                                    Triangle<>  ,
                                    Triangle<>  >::type ,
                      TetrahedronSideNodeMap >
{ typedef Tetrahedron<4> base ; };

template<> struct Tetrahedron<10> : public
  CellTopologyTraits< 3 , 4 , 10 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      TetrahedronEdgeNodeMap ,
                      MakeTypeList< Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  >::type ,
                      TetrahedronSideNodeMap >
{ typedef Tetrahedron<4> base ; };

template<> const CellTopology * cell_topology< Tetrahedron<> >();
template<> const CellTopology * cell_topology< Tetrahedron<10> >();

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

