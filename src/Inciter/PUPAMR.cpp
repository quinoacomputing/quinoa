// *****************************************************************************
/*!
  \file      src/Inciter/PUPAMR.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ Pack/UnPack utilities for AMR
  \details   This file contains some extensions to Charm++'s Pack/UnPack
    routines for use with AMR data structures.
*/
// *****************************************************************************

#include "PUPAMR.hpp"

void PUP::pup( PUP::er &p, AMR::Refinement_State& s )
// *****************************************************************************
//  Pack/Unpack Refinement_State
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] s Refinement_State object reference
// *****************************************************************************
{
  p | s.active_element_number;
  p | s.refinement_case;
  p | s.children;
  p | s.refinement_level;
  p | s.child_number;
  p | s.parent_id;
  p | s.normal;
}

void PUP::pup( PUP::er &p, AMR::Edge_Refinement& e )
// *****************************************************************************
//  Pack/Unpack Edge_Refinement
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] e Edge_Refinement object reference
// *****************************************************************************
{
  p | e.A;
  p | e.B;
  p | e.needs_refining;
  p | e.needs_derefining;
  p | e.lock_case;
}

void PUP::pup( PUP::er &p, AMR::edge_store_t& e )
// *****************************************************************************
//  Pack/Unpack edge_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] e edge_store_t object reference
// *****************************************************************************
{
  p | e.edges;
}

void PUP::pup( PUP::er &p, AMR::edge_t& e )
// *****************************************************************************
//  Pack/Unpack edge_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] e edge_t object reference
// *****************************************************************************
{
  p | e.get_data();
}

void PUP::pup( PUP::er &p, AMR::active_element_store_t& a )
// *****************************************************************************
//  Pack/Unpack active_element_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] a active_element_store_t object reference
// *****************************************************************************
{
  p | a.data();
}

void PUP::pup( PUP::er &p, AMR::master_element_store_t& m )
// *****************************************************************************
//  Pack/Unpack master_element_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] m master_element_store_t object reference
// *****************************************************************************
{
  p | m.data();
}

void PUP::pup( PUP::er &p, AMR::id_generator_t& i )
// *****************************************************************************
//  Pack/Unpack id_generator_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] i id_generator_t object reference
// *****************************************************************************
{
  p | i.start_id;
  p | i.next_tet_id;
}

void PUP::pup( PUP::er &p, AMR::tet_store_t& t )
// *****************************************************************************
//  Pack/Unpack tet_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] t tet_store_t object reference
// *****************************************************************************
{
  p | t.center_tets;
  p | t.delete_list;
  p | t.active_elements.data();
  p | t.master_elements.data();
  p | t.active_tetinpoel;
  p | t.active_nodes;
  p | t.id_generator;
  p | t.intermediate_list;
  p | t.active_id_mapping;
  p | t.tets;
  p | t.edge_store;
  p | t.marked_refinements;
  p | t.marked_derefinements;
}

void PUP::pup( PUP::er &p, AMR::mesh_adapter_t& m )
// *****************************************************************************
//  Pack/Unpack mesh_adapter_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] m mesh_adapter_t object reference
// *****************************************************************************
{
  p | m.derefinement_cut_off;
  p | m.refinement_cut_off;
  p | m.tet_store;
  p | m.node_connectivity;
#ifdef ENABLE_NODE_STORE
  p | m.node_store;
#endif
  p | m.refiner;
}

#ifdef ENABLE_NODE_STORE
void PUP::pup( PUP::er &p, AMR::node_store_t& n )
// *****************************************************************************
//  Pack/Unpack node_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] n node_store_t object reference
// *****************************************************************************
{
  p | n.m_x;
  p | n.m_y;
  p | n.m_z;
}
#endif


void PUP::pup( PUP::er &p, AMR::node_connectivity_t& n )
// *****************************************************************************
//  Pack/Unpack node_connectivity_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] n node_connectivity_t object reference
// *****************************************************************************
{
  p | n.data();
  p | n.inv_data();
  p | n.empty_node_count;
}

void PUP::pup( PUP::er &p, AMR::refinement_t& r )
// *****************************************************************************
//  Pack/Unpack refinement_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] r refinement_t object reference
// *****************************************************************************
{
  p | r.MAX_REFINEMENT_LEVEL;
}
