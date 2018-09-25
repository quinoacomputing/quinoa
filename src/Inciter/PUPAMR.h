// *****************************************************************************
/*!
  \file      src/Inciter/PUPAMR.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ Pack/UnPack utilities for AMR
  \details   This file contains some extensions to Charm++'s Pack/UnPack
    routines for use with AMR data structures.
*/
// *****************************************************************************
#ifndef PUPAMR_h
#define PUPAMR_h

#include "NoWarning/pup_stl.h"

#include "AMR/Refinement_State.h"
#include "AMR/edge_store.h"

//! Extensions to Charm++'s Pack/Unpack routines
namespace PUP {

//////////////////// Serialize Edge_Refinement ////////////////////

/** @name Charm++ pack/unpack serializer member functions */
///@{
//! Pack/Unpack Edge_Refinement
void pup( PUP::er &p, AMR::Edge_Refinement& e );

//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] e Edge_Refinement object reference
inline void operator|( PUP::er& p, AMR::Edge_Refinement& e ) { pup(p,e); }
//@}

//////////////////// Serialize edge_store_t ////////////////////

/** @name Charm++ pack/unpack serializer member functions */
///@{
//! Pack/Unpack edge_store_t
void pup( PUP::er &p, AMR::edge_store_t& e );

//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] e edge_store_t object reference
inline void operator|( PUP::er& p, AMR::edge_store_t& e ) { pup(p,e); }
//@}


} // PUP::

#endif // PUPAMR_h
