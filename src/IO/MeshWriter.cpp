// *****************************************************************************
/*!
  \file      src/IO/MeshWriter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ group for outputing mesh data to file
  \details   Charm++ group definition used to output data associated to
     unstructured meshes to file(s). Charm++ chares (work units) send mesh and
     field data associated to mesh entities to the MeshWriter class defined here
     to write the data to file(s).
*/
// *****************************************************************************

#include "QuinoaBuildConfig.hpp"
#include "MeshWriter.hpp"
#include "Reorder.hpp"
#include "ExodusIIMeshWriter.hpp"

using tk::MeshWriter;

MeshWriter::MeshWriter( ctr::FieldFileType filetype,
                        Centering bnd_centering,
                        bool benchmark,
                        std::size_t nmesh ) :
  m_filetype( filetype ),
  m_bndCentering( bnd_centering ),
  m_benchmark( benchmark ),
  m_nmesh( nmesh ),
  m_nchare( nmesh, 0 )
// *****************************************************************************
//  Constructor: set some defaults that stay constant at all times
//! \param[in] filetype Output file format type
//! \param[in] bnd_centering Centering to identify what boundary data to write.
//!   For a nodal scheme, e.g., ALECG, this is nodal, for a DG scheme, this is
//!   cell-based.
//! \param[in] benchmark True of benchmark mode. No field output happens in
//!   benchmark mode. This (and associated if tests) are here so client code
//!   does not have to deal with this.
//! \param[in] nmesh Total number of meshes
// *****************************************************************************
{
}

void
MeshWriter::nchare( std::size_t meshid, int n )
// *****************************************************************************
//  Set the total number of chares
//! \param[in] meshid Mesh whose number of chares to set
//! \param[in] n Total number of chares across the whole problem (for a mesh)
// *****************************************************************************
{
  Assert( meshid < m_nchare.size(), "Size mismatch" );
  m_nchare[ meshid ] = n;
}

void
MeshWriter::write(
  std::size_t meshid,
  bool meshoutput,
  bool fieldoutput,
  uint64_t itr,
  uint64_t itf,
  tk::real time,
  int chareid,
  const std::string& basefilename,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel,
  const std::vector< std::string >& elemfieldnames,
  const std::vector< std::string >& nodefieldnames,
  const std::vector< std::string >& elemsurfnames,
  const std::vector< std::string >& nodesurfnames,
  const std::vector< std::vector< tk::real > >& elemfields,
  const std::vector< std::vector< tk::real > >& nodefields,
  const std::vector< std::vector< tk::real > >& elemsurfs,
  const std::vector< std::vector< tk::real > >& nodesurfs,
  const std::set< int >& outsets,
  CkCallback c )
// *****************************************************************************
//  Output unstructured mesh into file
//! \param[in] meshid Mesh Id
//! \param[in] meshoutput True if mesh is to be written
//! \param[in] fieldoutput True if field data is to be written
//! \param[in] itr Iteration count since a new mesh. New mesh in this context
//!   means that either the mesh is moved and/or its topology has changed.
//! \param[in] itf Field output iteration count
//! \param[in] time Physical time this at this field output dump
//! \param[in] chareid The chare id the write-to-file request is coming from
//! \param[in] basefilename String to use as the base of the filename
//! \param[in] inpoel Mesh connectivity for the mesh chunk to be written with
//!   local ids
//! \param[in] coord Node coordinates of the mesh chunk to be written
//! \param[in] bface Map of boundary-face lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] bnode Map of boundary-node lists mapped to corresponding side set
//!   ids for this mesh chunk with local ids
//! \param[in] triinpoel Interconnectivity of points and boundary-face in this
//!   mesh chunk with local ids
//! \param[in] elemfieldnames Names of element fields to be output to file
//! \param[in] nodefieldnames Names of node fields to be output to file
//! \param[in] elemsurfnames Names of elemental surface fields to be output to
//!   file
//! \param[in] nodesurfnames Names of node surface fields to be output to file
//! \param[in] elemfields Field data in mesh elements to output to file
//! \param[in] nodefields Field data in mesh nodes to output to file
//! \param[in] elemsurfs Surface field data in mesh elements to output to file
//! \param[in] nodesurfs Surface field data in mesh nodes to output to file
//! \param[in] outsets Unique set of surface side set ids along which to save
//!   solution field variables
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  if (!m_benchmark) {

    // Generate filenames for volume and surface field output
    auto vf = filename( basefilename, meshid, itr, chareid );
  
    if (meshoutput) {
      if (m_filetype == ctr::FieldFileType::EXODUSII) {

        // Write volume mesh and field names
        ExodusIIMeshWriter ev( vf, ExoWriter::CREATE );
        // Write chare mesh (do not write side sets in parallel)
        if (m_nchare[meshid] == 1) {

          if (m_bndCentering == Centering::ELEM)
            ev.writeMesh( inpoel, coord, bface, triinpoel );
          else if (m_bndCentering == Centering::NODE)
            ev.writeMesh( inpoel, coord, bnode );
          else Throw( "Centering not handled for writing mesh" );

        } else {
          ev.writeMesh< 4 >( inpoel, coord );
        }
        ev.writeElemVarNames( elemfieldnames );
        Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );
        ev.writeNodeVarNames( nodefieldnames );

        // Write surface meshes and surface variable field names
        for (auto s : outsets) {
          auto sf = filename( basefilename, meshid, itr, chareid, s );
          ExodusIIMeshWriter es( sf, ExoWriter::CREATE );
          auto b = bface.find(s);
          if (b == end(bface)) {
            // If a side set does not exist on a chare, write out a
            // connectivity for a single triangle with its node coordinates of
            // zero. This is so the paraview series reader can load side sets
            // distributed across multiple files. See also
            // https://www.paraview.org/Wiki/Restarted_Simulation_Readers.
            es.writeMesh< 3 >( std::vector< std::size_t >{1,2,3},
              UnsMesh::Coords{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }} );
            es.writeElemVarNames( elemsurfnames );
            es.writeNodeVarNames( nodesurfnames );
            continue;
          }
          std::vector< std::size_t > nodes;
          for (auto f : b->second) {
            nodes.push_back( triinpoel[f*3+0] );
            nodes.push_back( triinpoel[f*3+1] );
            nodes.push_back( triinpoel[f*3+2] );
          }
          auto [inp,gid,lid] = tk::global2local( nodes );
          tk::unique( nodes );
          auto nnode = nodes.size();
          UnsMesh::Coords scoord;
          scoord[0].resize( nnode );
          scoord[1].resize( nnode );
          scoord[2].resize( nnode );
          std::size_t j = 0;
          for (auto i : nodes) {
            scoord[0][j] = coord[0][i];
            scoord[1][j] = coord[1][i];
            scoord[2][j] = coord[2][i];
            ++j;
          }
          es.writeMesh< 3 >( inp, scoord );
          es.writeElemVarNames( elemsurfnames );
          es.writeNodeVarNames( nodesurfnames );
        }

      }
    }

    if (fieldoutput) {
      if (m_filetype == ctr::FieldFileType::EXODUSII) {

        // Write volume variable fields
        ExodusIIMeshWriter ev( vf, ExoWriter::OPEN );
        ev.writeTimeStamp( itf, time );
        // Write volume element variable fields
        int varid = 0;
        for (const auto& v : elemfields) ev.writeElemScalar( itf, ++varid, v );
        // Write volume node variable fields
        varid = 0;
        for (const auto& v : nodefields) ev.writeNodeScalar( itf, ++varid, v );

        // Write surface node variable fields
        std::size_t j = 0;
        std::size_t k = 0;
        auto nvar = static_cast< int >( nodesurfnames.size() ) ;
        auto nevar = static_cast< int >( elemsurfnames.size() ) ;
        for (auto s : outsets) {
          auto sf = filename( basefilename, meshid, itr, chareid, s );
          ExodusIIMeshWriter es( sf, ExoWriter::OPEN );
          es.writeTimeStamp( itf, time );
          if (bface.find(s) == end(bface)) {
            // If a side set does not exist on a chare, write out a
            // a node field for a single triangle with zeros. This is so the
            // paraview series reader can load side sets distributed across
            // multiple files. See also
            // https://www.paraview.org/Wiki/Restarted_Simulation_Readers.
            for (int i=1; i<=nvar; ++i) es.writeNodeScalar( itf, i, {0,0,0} );
            for (int i=1; i<=nevar; ++i) es.writeElemScalar( itf, i, {0} );
            continue;
          }
          for (int i=1; i<=nvar; ++i)
            es.writeNodeScalar( itf, i, nodesurfs[j++] );
          for (int i=1; i<=nevar; ++i)
            es.writeElemScalar( itf, i, elemsurfs[k++] );
        }

      }
    }

  }

  c.send();
}

std::string
MeshWriter::filename( const std::string& basefilename,
                      std::size_t meshid,
                      uint64_t itr,
                      int chareid,
                      int surfid ) const
// *****************************************************************************
//  Compute filename
//! \param[in] basefilename String use as the base filename.
//! \param[in] meshid Mesh Id
//! \param[in] itr Iteration count since a new mesh. New mesh in this context
//!   means that either the mesh is moved and/or its topology has changed.
//! \param[in] chareid The chare id the write-to-file request is coming from
//! \param[in] surfid Surface ID if computing a surface filename
//! \details We use a file naming convention for large field output data that
//!   allows ParaView to glue multiple files into a single simulation output by
//!   only loading a single file. The base filename is followed by ".e-s.",
//!   which probably stands for Exodus Sequence, followed by 3 integers:
//!   (1) {RS}: counts the number of "restart dumps", but we use this for
//!   counting the number of outputs with a different mesh, e.g., due to
//!   mesh refinement, thus if this first number is new the mesh is new
//!   compared to the previous (first) number afer ".e-s.",
//!   (2) {NP}: total number of partitions (workers, chares), this is more than
//!   the number of PEs with nonzero virtualization (overdecomposition), and
//!   (3) {RANK}: worker (chare) id.
//!   Thus {RANK} does spatial partitioning, while {RS} partitions in time, but
//!   a single {RS} id may contain multiple time steps, which equals to the
//!   number of time steps at which field output is saved without refining the
//!   mesh.
//! \return Filename computed
//! \see https://www.paraview.org/Wiki/Restarted_Simulation_Readers
// *****************************************************************************
{
  return basefilename + (surfid ? "-surf." + std::to_string(surfid) : "")
         + (m_nmesh > 1 ? '.' + std::to_string(meshid) : "")
         + ".e-s"
         + '.' + std::to_string( itr )        // iteration count with new mesh
         + '.' + std::to_string( m_nchare[meshid] ) // total number of workers
         + '.' + std::to_string( chareid )    // new file per worker
         ;
}

#include "NoWarning/meshwriter.def.h"
