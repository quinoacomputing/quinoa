// *****************************************************************************
/*!
  \file      src/Inciter/PUPAMR.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ Pack/UnPack utilities for AMR
  \details   This file contains some extensions to Charm++'s Pack/UnPack
    routines for use with AMR data structures.
*/
// *****************************************************************************

#include "PUPAMR.h"

void PUP::pup( PUP::er &p, AMR::Edge_Refinement& e )
// *****************************************************************************
//  Pack/Unpack Edge_Refinement
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] e Edge_Refinement object reference
// *****************************************************************************
{
  p | e.A;
  p | e.B;
  p | e.refinement_criteria;
  p | e.needs_refining;
  p | e.needs_derefining;
  p | e.is_dead;
  p | e.lock_case;
}

void PUP::pup( PUP::er &p, AMR::edge_store_t& e )
// *****************************************************************************
//  Pack/Unpack edge_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] e edge_store_T object reference
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

void PUP::pup( PUP::er &p, AMR::marked_refinements_store_t& m )
// *****************************************************************************
//  Pack/Unpack marked_refinements_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] m marked_refinements_store_t object reference
// *****************************************************************************
{
  p | m.data();
}

void PUP::pup( PUP::er &p, AMR::tet_store_t& t )
// *****************************************************************************
//  Pack/Unpack tet_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] t tet_store_t object reference
// *****************************************************************************
{
  p | t.intermediate_list;
  p | t.active_id_mapping;
  p | t.tets;
  p | t.marked_refinements;
}

void PUP::pup( PUP::er &p, AMR::mesh_adapter_t& m )
// *****************************************************************************
//  Pack/Unpack mesh_adapter_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] m mesh_adapter_t object reference
// *****************************************************************************
{
  p | m.tet_store;
  p | m.node_connectivity;
  p | m.node_store;
  p | m.refiner;
}

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
  p | n.m_graphsize;
}

void PUP::pup( PUP::er &p, AMR::node_connectivity_t& n )
// *****************************************************************************
//  Pack/Unpack node_connectivity_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] n node_connectivity_t object reference
// *****************************************************************************
{
  p | n.data();
}

void PUP::pup( PUP::er &, AMR::refinement_t& )
// *****************************************************************************
//  Pack/Unpack refinement_t
// *****************************************************************************
{
  // no state in refinement_t currently
}
