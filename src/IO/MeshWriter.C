// *****************************************************************************
/*!
  \file      src/IO/MeshWriter.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ group for outputing mesh data to file
  \details   Charm++ group definition used to output data associated to
     unstructured meshes to file(s). Charm++ chares (work units) send mesh and
     field data associated to mesh entities to the MeshWriter class defined here
     to write the data to file(s).
*/
// *****************************************************************************

#include "QuinoaConfig.h"
#include "MeshWriter.h"
#include "ExodusIIMeshWriter.h"

#ifdef HAS_ROOT
  #include "RootMeshWriter.h"
#endif

using tk::MeshWriter;

MeshWriter::MeshWriter( ctr::FieldFileType filetype,
                        Centering bnd_centering,
                        bool benchmark ) :
  m_filetype( filetype ),
  m_bndCentering( bnd_centering ),
  m_benchmark( benchmark ),
  m_nchare( 0 )
// *****************************************************************************
//  Constructor: set some defaults that stay constant at all times
//! \param[in] filetype Output file format type
//! \param[in] bnd_centering Centering to identify what boundary data to write.
//!   For a nodal scheme, e.g., DiagCG, this is nodal, for a DG scheme, this is
//!   cell-based.
//! \param[in] benchmark True of benchmark mode. No field output happens in
//!   benchmark mode. This (and associated if tests) are here so client code
//!   does not have to deal with this.
// *****************************************************************************
{
}

void
MeshWriter::nchare( int n )
// *****************************************************************************
//  Set the total number of chares
//! \param[in] n Total number of chares across the whole problem
// *****************************************************************************
{
  m_nchare = n;
}

void
MeshWriter::write(
  bool meshoutput,
  bool fieldoutput,
  uint64_t itr,
  uint64_t itf,
  tk::real time,
  int chareid,
  Centering centering,
  const std::string& basefilename,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const std::vector< std::string >& names,
  const std::vector< std::vector< tk::real > >& fields,
  CkCallback c )
// *****************************************************************************
//  Output unstructured mesh into file
//! \param[in] meshoutput True if mesh is to be written
//! \param[in] fieldoutput True if field data is to be written
//! \param[in] itr Iteration count since a new mesh. New mesh in this context
//!   means that either the mesh is moved and/or its topology has changed.
//! \param[in] itf Field output iteration count
//! \param[in] time Physical time this at this field output dump
//! \param[in] chareid The chare id the write-to-file request is coming from
//! \param[in] basefilename String to use as the base of the filename
//! \param[in] centering The centering that will be associated to the field data
//!   to be output when writeFields is called next
//! \param[in] inpoel Mesh connectivity for the mesh chunk to be written
//! \param[in] coord Node coordinates of the mesh chunk to be written
//! \param[in] bface Map of boundary-face lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] triinpoel Interconnectivity of points and boundary-face in this
//!   mesh chunk
//! \param[in] bnode Map of boundary-node lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] lid Global->local node id map for the mesh chunk to be written
//! \param[in] names Names of fields to be output in next call to writeFields()
//! \param[in] fields Field data to output to file
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  if (!m_benchmark) {

    auto f = filename( basefilename, itr, chareid );
  
    if (meshoutput) {
      #ifdef HAS_ROOT
      if (m_filetype == ctr::FieldFileType::ROOT) {

        RootMeshWriter rmw( f, 0 );
        rmw.writeMesh( UnsMesh( inpoel, coord ) );
        rmw.writeNodeVarNames( names );

      } else
      #endif
      if (m_filetype == ctr::FieldFileType::EXODUSII) {
        ExodusIIMeshWriter ew( f, ExoWriter::CREATE );
        // Write chare mesh (do not write side sets in parallel)
        if (m_nchare == 1) {

          if (m_bndCentering == Centering::ELEM)
            ew.writeMesh( inpoel, coord, bface, triinpoel );
          else if (m_bndCentering == Centering::NODE) {
            // Convert boundary node lists to local ids for output
            std::map< int, std::vector< std::size_t > > lbnode = bnode;
            for (auto& s : lbnode)
              for (auto& p : s.second)
                p = cref_find( lid, p );
            ew.writeMesh( inpoel, coord, lbnode );
          } else Throw( "Centering not handled for writing mesh" );

        } else {
          ew.writeMesh( inpoel, coord );

        }
        // Write field names
        if (centering == Centering::ELEM)
          ew.writeElemVarNames( names );
        else if (centering == Centering::NODE)
          ew.writeNodeVarNames( names );
      }
    }

    if (fieldoutput) {
      #ifdef HAS_ROOT
      if (m_filetype == ctr::FieldFileType::ROOT) {

        RootMeshWriter rw( f, 1 );
        rw.writeTimeStamp( itf, time );
        int varid = 0;
        for (const auto& v : fields) rw.writeNodeScalar( itf, ++varid, v );

      } else
      #endif
      if (m_filetype == ctr::FieldFileType::EXODUSII) {

        ExodusIIMeshWriter ew( f, ExoWriter::OPEN );
        ew.writeTimeStamp( itf, time );
        if (centering == Centering::ELEM) {
          int varid = 0;
          for (const auto& v : fields) ew.writeElemScalar( itf, ++varid, v );
        } else if (centering == Centering::NODE) {
          int varid = 0;
          for (const auto& v : fields) ew.writeNodeScalar( itf, ++varid, v );
        }

      }
    }

  }

  c.send();
}

std::string
MeshWriter::filename( const std::string& basefilename,
                      uint64_t itr,
                      int chareid ) const
// *****************************************************************************
//  Compute filename
//! \param[in] itr Iteration count since a new mesh. New mesh in this context
//!   means that either the mesh is moved and/or its topology has changed.
//! \param[in] chareid The chare id the write-to-file request is coming from
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
  return basefilename + ".e-s"
         + '.' + std::to_string( itr )        // iteration count with new mesh
         + '.' + std::to_string( m_nchare )   // total number of workers
         + '.' + std::to_string( chareid )    // new file per worker
         #ifdef HAS_ROOT
         + (m_filetype == ctr::FieldFileType::ROOT ? ".root" : "")
         #endif
         ;
}

#include "NoWarning/meshwriter.def.h"
