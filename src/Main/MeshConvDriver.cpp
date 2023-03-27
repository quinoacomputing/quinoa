// *****************************************************************************
/*!
  \file      src/Main/MeshConvDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
// *****************************************************************************

#include "Types.hpp"
#include "Tags.hpp"
#include "MeshConvDriver.hpp"
#include "MeshFactory.hpp"
#include "TaggedTupleDeepPrint.hpp"
#include "Writer.hpp"

#include "NoWarning/meshconv.decl.h"

using meshconv::MeshConvDriver;

extern CProxy_Main mainProxy;

MeshConvDriver::MeshConvDriver( const ctr::CmdLine& cmdline, int ) :
  m_print( cmdline.logname( cmdline.get< tag::io, tag::screen >(),
                            cmdline.get< tag::io, tag::nrestart >() ),
           cmdline.get< tag::verbose >() ? std::cout : std::clog,
           std::ios_base::app ),
  m_reorder( cmdline.get< tag::reorder >() ),
  m_input(),
  m_output()
// *****************************************************************************
//  Constructor
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
// *****************************************************************************
{
  // Save input file name
  m_input = cmdline.get< tag::io, tag::input >();
  // Save output file name
  m_output = cmdline.get< tag::io, tag::output >();

  // Output command line object to file
  auto logfilename = tk::meshconv_executable() + "_input.log";
  tk::Writer log( logfilename );
  tk::print( log.stream(), "cmdline", cmdline );
}

void
MeshConvDriver::execute() const
// *****************************************************************************
//  Execute: Convert mesh file
// *****************************************************************************
{
  m_print.endsubsection();

  std::vector< std::pair< std::string, tk::real > > times;

  // If input filename contains a '%', we aggregate multiple files
  if (m_input.find('%') == std::string::npos) {

    // Convert single mesh

    times.push_back( {} );
    auto mesh = tk::readUnsMesh( m_print, m_input, times[0] );
    auto wtimes = tk::writeUnsMesh( m_print, m_output, mesh, m_reorder );
    times.insert( end(times), begin(wtimes), end(wtimes) );

  } else {

    // Aggregate multiple meshes containing surface output

    // Find a '%' sign in the input filename, and assuming a syntax of
    // '.<nfile>.%', find '<nfile>' as the number of files to aggregate.
    auto percent_pos = m_input.find( '%' );
    auto input_basename = m_input.substr( 0, percent_pos );
    auto dot1 = m_input.find_last_of( '.', percent_pos );
    auto dot2 = m_input.find_last_of( '.', dot1-1 );
    auto nfile_str = m_input.substr( dot2+1, dot1-dot2-1  );
    std::stringstream ss( nfile_str );
    if (nfile_str.empty())
      Throw( "The percent sign must be followed by an "
             "integer, the number of files to aggregate" );
    std::size_t nfile;
    ss >> nfile;
    m_print.diag( "Aggregating " + std::to_string(nfile) +
                  " files from base filename: '" + input_basename +'\'' );

    const auto eps = std::numeric_limits< tk::real >::epsilon();

    // Lambda to echo some diagnostics on the mesh being processes to screen
    auto diag = [&]( const std::string& name, const tk::UnsMesh& mesh ){
      m_print.diag( name + ": ntri: " +
        std::to_string(mesh.triinpoel().size()/3) +
        ", ntime: " + std::to_string(mesh.vartimes().size()) +
        (!mesh.nodevars().empty() ? ", node_var: " +
           std::to_string(mesh.nodevars()[0].size()) : "") +
        (!mesh.nodevars()[0].empty() ? ", npoin: " +
           std::to_string(mesh.nodevars()[0][0].size()) : "") +
        (!mesh.elemvars().empty() ? ", elem_var: " +
           std::to_string(mesh.elemvars()[0].size()) : "") +
        (!mesh.elemvars()[0].empty() ? ", nelem: " +
           std::to_string(mesh.elemvars()[0][0].size()) : "") );
    };

    // Output-mesh containers, will store aggregated surface(s) and field output
    tk::UnsMesh::Coords coords;
    auto& X = coords[0];
    auto& Y = coords[1];
    auto& Z = coords[2];
    std::size_t npoin = 0;
    std::size_t nelem = 0;
    std::vector< std::size_t > otriinpoel;
    std::vector< std::string > nodevarnames;
    std::vector< std::string > elemvarnames;
    std::vector< tk::real > vartimes;
    std::vector< std::vector< std::vector< tk::real > > > nodevars;
    std::vector< std::vector< std::vector< tk::real > > > elemvars;
    // Counter for number of non-empty meshes processed
    std::size_t k = 0;
    for (std::size_t m=0; m<nfile; ++m) {
      std::string name = input_basename + std::to_string(m);
      times.push_back( {} );
      auto mesh = tk::readUnsMesh( m_print, name, times.back() );
      const auto& triinpoel = mesh.triinpoel();
      // Skip meshes with a single triange cell
      if (triinpoel.size() == 3) continue;
      const auto& x = mesh.x();
      const auto& y = mesh.y();
      const auto& z = mesh.z();
      nodevarnames = mesh.nodevarnames();
      elemvarnames = mesh.elemvarnames();
      vartimes = mesh.vartimes();
      // Echo some diagnostics on the mesh being processes to screen
      diag( name, mesh );
      // Aggregate data from each triangle element in mesh
      for (std::size_t e=0; e<triinpoel.size()/3; ++e) {
        for (std::size_t n=0; n<3; ++n) {
          auto j = triinpoel[ e*3+n ];
          bool visited = false;
          // WARNING: linear search below, will not scale well
          for (std::size_t i=0; i<X.size(); ++i) {
            // If mesh point has already been seen (on a previous mesh)
            if (std::abs(x[j]-X[i]) < eps &&
                std::abs(y[j]-Y[i]) < eps &&
                std::abs(z[j]-Z[i]) < eps)
            { // no point in connectivity but nothing else
              visited = true;
              otriinpoel.push_back( i );
            }
          }
          if (!visited) { // Mesh point not yet seen
            // save coordinates and (global) point id in aggregated connectivity
            X.push_back( x[j] );
            Y.push_back( y[j] );
            Z.push_back( z[j] );
            otriinpoel.push_back( npoin );
            // aggregate nodal field data for all times and variables
            std::size_t time = 0;
            std::size_t varid = 0;
            for (const auto& t : mesh.nodevars()) {  // for all times
              if (k == 0 && npoin == 0) nodevars.push_back( {} );
              for (const auto& v : t) {              // for all variables
                if (k == 0 && npoin == 0) nodevars.back().push_back( {} );
                nodevars[time][varid].push_back( v[j] );
                ++varid;
              }
              ++time;
              varid = 0;
            }
            ++npoin;      // increase number of nodes in output mesh
          }
        }

        // aggregate elemental field data for all times and variables
        std::size_t etime = 0;
        std::size_t evarid = 0;
        for (const auto& t : mesh.elemvars()) {  // for all times
          if (k == 0 && nelem == 0) elemvars.push_back( {} );
          for (const auto& v : t) {              // for all variables
            if (k == 0 && nelem == 0) elemvars.back().push_back( {} );
            elemvars[etime][evarid].push_back( v[e] );
            ++evarid;
          }
          ++etime;
          evarid = 0;
        }
        ++nelem;    // increase number of elements in output mesh
      }
      ++k;        // increase number of non-empty meshes processed
    }

    // Construct aggregated output mesh
    tk::UnsMesh outmesh( coords, otriinpoel, nodevarnames, elemvarnames,
      vartimes, nodevars, elemvars );
    // Echo diagnostics on the aggreegate output mesh
    diag( m_output, outmesh );
    // Write output mesh to file
    auto wtimes = tk::writeUnsMesh( m_print, m_output, outmesh, m_reorder );
    // Collect wall-clock time data
    times.insert( end(times), begin(wtimes), end(wtimes) );

  }

  mainProxy.timestamp( times );

  mainProxy.finalize();
}
