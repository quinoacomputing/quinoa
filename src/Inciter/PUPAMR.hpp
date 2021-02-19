// *****************************************************************************
/*!
  \file      src/Inciter/PUPAMR.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ Pack/UnPack utilities for AMR
  \details   This file contains some extensions to Charm++'s Pack/UnPack
    routines for use with AMR data structures.
*/
// *****************************************************************************
#ifndef PUPAMR_h
#define PUPAMR_h

#include "Base/PUPUtil.hpp"

#include "AMR/edge_store.hpp"
#include "AMR/edge.hpp"
#include "AMR/marked_refinements_store.hpp"
#include "AMR/tet_store.hpp"
#include "AMR/mesh_adapter.hpp"
#include "AMR/node_store.hpp"
#include "AMR/node_connectivity.hpp"
#include "AMR/refinement.hpp"
#include "AMR/master_element_store.hpp"
#include "AMR/id_generator.hpp"

//! Extensions to Charm++'s Pack/Unpack routines
namespace PUP {

/** @name Charm++ pack/unpack serializer member functions for Refinement_State */
///@{
//! Pack/Unpack Refinement_State
void pup( PUP::er &p, AMR::Refinement_State& s );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] s Refinement_State object reference
inline void operator|( PUP::er& p, AMR::Refinement_State& s ) { pup(p,s); }
//@}

/** @name Charm++ pack/unpack serializer member functions for Edge_Refinement */
///@{
//! Pack/Unpack Edge_Refinement
void pup( PUP::er &p, AMR::Edge_Refinement& e );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] e Edge_Refinement object reference
inline void operator|( PUP::er& p, AMR::Edge_Refinement& e ) { pup(p,e); }
//@}

/** @name Charm++ pack/unpack serializer member functions for edge_store_t */
///@{
//! Pack/Unpack edge_store_t
void pup( PUP::er &p, AMR::edge_store_t& e );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] e edge_store_t object reference
inline void operator|( PUP::er& p, AMR::edge_store_t& e ) { pup(p,e); }
//@}

/** @name Charm++ pack/unpack serializer member functions for edge_t */
///@{
//! Pack/Unpack edge_t
void pup( PUP::er &p, AMR::edge_t& e );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] e edge_t object reference
inline void operator|( PUP::er& p, AMR::edge_t& e ) { pup(p,e); }
//@}

/** @name Charm++ pack/unpack serializer member functions for marked_refinements_store_t */
///@{
//! Pack/Unpack marked_refinements_store_t
//! \param[in] p Charm++'s pack/unpack object
//! \param[in,out] m marked_refinements_store_t object reference
template< class case_t >
void pup( PUP::er &p, AMR::marked_refinements_store_t< case_t >& m ) {
  p | m.data();
  p | m.get_state_changed();
}
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] m marked_refinements_store_t object reference
template< class case_t >
inline void operator|( PUP::er& p, AMR::marked_refinements_store_t<case_t>& m )
{ pup(p,m); }
//@}

/** @name Charm++ pack/unpack serializer member functions for active_element_store_t */
///@{
//! Pack/Unpack active_element_store_t
void pup( PUP::er &p, AMR::active_element_store_t& a );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] a active_element_store_t object reference
inline void operator|( PUP::er& p, AMR::active_element_store_t& a )
{ pup(p,a); }
//@}

/** @name Charm++ pack/unpack serializer member functions for master_element_store_t */
///@{
//! Pack/Unpack master_element_store_t
void pup( PUP::er &p, AMR::master_element_store_t& m );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] m master_element_store_t object reference
inline void operator|( PUP::er& p, AMR::master_element_store_t& m )
{ pup(p,m); }
//@}
/** @name Charm++ pack/unpack serializer member functions for active_element_store_t */

///@{
//! Pack/Unpack id_generator_t
void pup( PUP::er &p, AMR::id_generator_t& i );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] i id_generator_t object reference
inline void operator|( PUP::er& p, AMR::id_generator_t& i )
{ pup(p,i); }
//@}

/** @name Charm++ pack/unpack serializer member functions for tet_store_t */
///@{
//! Pack/Unpack tet_store_t
void pup( PUP::er &p, AMR::tet_store_t& t );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] t tet_store_t object reference
inline void operator|( PUP::er& p, AMR::tet_store_t& t ) { pup(p,t); }
//@}

/** @name Charm++ pack/unpack serializer member functions for mesh_adapter_t */
///@{
//! Pack/Unpack mesh_adapter_t
void pup( PUP::er &p, AMR::mesh_adapter_t& m );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] m mesh_adapter_t object reference
inline void operator|( PUP::er& p, AMR::mesh_adapter_t& m ) { pup(p,m); }
//@}

/** @name Charm++ pack/unpack serializer member functions for node_store_t */
///@{
//! Pack/Unpack node_store_t
void pup( PUP::er &p, AMR::node_store_t& n );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] n node_store_t object reference
inline void operator|( PUP::er& p, AMR::node_store_t& n ) { pup(p,n); }
//@}

/** @name Charm++ pack/unpack serializer member functions for node_connectivity_t */
///@{
//! Pack/Unpack node_connectivity_t
void pup( PUP::er &p, AMR::node_connectivity_t& n );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] n node_connectivity_t object reference
inline void operator|( PUP::er& p, AMR::node_connectivity_t& n ) { pup(p,n); }
//@}

/** @name Charm++ pack/unpack serializer member functions for refinement_t */
///@{
//! Pack/Unpack refinement_t
void pup( PUP::er &p, AMR::refinement_t& r );
//! Pack/Unpack serialize operator|
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] r refinement_t object reference
inline void operator|( PUP::er& p, AMR::refinement_t& r ) { pup(p,r); }
//@}

} // PUP::

#endif // PUPAMR_h
