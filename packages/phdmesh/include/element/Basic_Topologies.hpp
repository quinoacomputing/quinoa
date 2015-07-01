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

#ifndef phdmesh_Basic_Topologies_hpp
#define phdmesh_Basic_Topologies_hpp

#include <element/CellTopology.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
// Node   < Dim = 0, Vertices = 0, Nodes = 0 >
// Line   < Dim = 1, Vertices = 2, Nodes = 2 >
//
// Line_3 < Dim = 1, Vertices = 2, Nodes = 3 >
//
//     [0]----[3]----[2]----[4]----[1]  ---> Positive direction
//----------------------------------------------------------------------

typedef CellTopologyTraits<0,0,0> Node_Traits ;

template< unsigned Nodes = 2 > struct Line ;
template< unsigned Nodes = 2 > struct ShellLine ;

template<> struct Line<2> : public CellTopologyTraits<1,2,2>
{ typedef Line<2> base ; };

template<> struct Line<3> : public CellTopologyTraits<1,2,3>
{ typedef Line<2> base ; };


typedef
  MakeTypeList< IndexList< 0 , 1 , 2 > ,
                IndexList< 1 , 0 , 2 > >::type
    ShellLineEdgeNodeMap ;

template<> struct ShellLine<2> : public
  CellTopologyTraits< 2 , 2 , 2 ,
                      MakeTypeList< Line<> , Line<> >::type ,
                      ShellLineEdgeNodeMap >
{ typedef ShellLine<2> base ; };

template<> struct ShellLine<3> : public
  CellTopologyTraits< 2 , 2 , 3 ,
                      MakeTypeList< Line<> , Line<> >::type ,
                      ShellLineEdgeNodeMap >
{ typedef ShellLine<2> base ; };

template<> const CellTopology * cell_topology< Node_Traits >();

template<> const CellTopology * cell_topology< Line<2> >();
template<> const CellTopology * cell_topology< Line<3> >();

template<> const CellTopology * cell_topology< ShellLine<2> >();
template<> const CellTopology * cell_topology< ShellLine<3> >();

} // namespace phdmesh

#endif

