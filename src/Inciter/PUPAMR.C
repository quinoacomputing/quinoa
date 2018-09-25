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

#include "NoWarning/charm++.h"

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
