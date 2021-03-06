// *****************************************************************************
/*!
  \file      src/PDE/InfoMesh.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Assemble info on the configuration of solver coupling
  \details   Assemble info on the configuration of solver coupling.
*/
// *****************************************************************************
#ifndef InfoMesht_h
#define InfoMesht_h

#include "ContainerUtil.hpp"

namespace inciter {

//! Extract configuration info on mesh (solver) coupling
//! \tparam eq Equation (solver) type
//! \param[in] c Index of eq (solver) to extract info on
//! \param[in,out] nfo Info object to augment
//! \note The index c is the index of the solver within the solvers of type eq.
template< class eq >
void infoMesh( std::size_t c,
               std::vector< std::pair< std::string, std::string > >& nfo )
{
  using tk::parameters;

  const auto& mesh = g_inputdeck.template get< tag::param, eq, tag::mesh >();
  const auto& mesh_filename = mesh.template get< tag::filename >();

  if (mesh_filename.size() > c) {
    nfo.emplace_back( "mesh id",
      std::to_string( mesh.template get< tag::id >()[c] ) );
    nfo.emplace_back( "mesh", mesh.template get< tag::filename >()[c] );
  }

  const auto& mesh_reference = mesh.template get< tag::reference >();

  if (mesh_reference.size() > c && mesh_reference[c] != '-') {
    nfo.emplace_back( "mesh reference", std::string( 1, mesh_reference[c] ) );
    const auto& mesh_location = mesh.template get< tag::location >();
    if (mesh_location.size() > c)
      nfo.emplace_back( "mesh location", parameters( mesh_location[c] ) );
    const auto& mesh_orientation = mesh.template get< tag::orientation >();
    if (mesh_orientation.size() > c)
      nfo.emplace_back( "mesh orientation", parameters( mesh_orientation[c] ) );
  }
}

} // inciter::

#endif // InfoMesh_h
