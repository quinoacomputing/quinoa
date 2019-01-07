// *****************************************************************************
/*!
  \file      src/IO/MeshWriter.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ group for outputing mesh data to file
  \details   Charm++ group definition used to output data associated to
     unstructured meshes to file(s). Charm++ chares (work units) send mesh and
     field data associated to mesh entities to the MeshWriter class defined here
     to write the data to file(s).
*/
// *****************************************************************************

#include <iostream>     // NOT NEEDED

#include "QuinoaConfig.h"
#include "MeshWriter.h"
#include "ExodusIIMeshWriter.h"

#ifdef HAS_ROOT
  #include "RootMeshWriter.h"
#endif

using tk::MeshWriter;

MeshWriter::MeshWriter( const std::string& output_basefilename,
                        ctr::FieldFileType filetype,
                        Centering bnd_centering,
                        bool benchmark ) :
  m_outputBasefilename( output_basefilename ),
  m_filetype( filetype ),
  m_bndCentering( bnd_centering ),
  m_benchmark( benchmark ),
  m_nchare( 0 )
// *****************************************************************************
//  Constructor: set some defaults that stay constant at all times
//! \param[in] output_basefilename String to use as the base of the filename
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
MeshWriter::expect( int nchare )
// *****************************************************************************
//  Set the total number of chares
//! \param[in] nchare Total number of chares across the whole problem
// *****************************************************************************
{
  m_nchare = nchare;
}

void
MeshWriter::writeMesh(
  uint64_t itr,
  int chareid,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::unordered_map< std::size_t, std::size_t >& lid )
// *****************************************************************************
//  Output unstructured mesh into file
//! \param[in] itr Iteration count since a new mesh (New mesh in this context
//!   means, either the mesh is moved and/or its topology has changed
//! \param[in] chareid The chare id the write-to-file request is coming from
//! \param[in] inpoel Mesh connectivity for the mesh chunk to be written
//! \param[in] coord Node coordinates of the mesh chunk to be written
//! \param[in] bface Map of boundary-face lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] triinpoel Interconnectivity of points and boundary-face in this
//!   mesh chunk
//! \param[in] bnode Map of boundary-node lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] lid Global->local node id map for the mesh chunk to be written
// *****************************************************************************
{
  if (!m_benchmark) {

    auto f = filename( itr, chareid );
  
    #ifdef HAS_ROOT
    if (m_filetype == ctr::FieldFileType::ROOT) {
  
      RootMeshWriter rmw( f, 0 );
      rmw.writeMesh( UnsMesh( inpoel, coord ) );
  
    } else
    #endif
    if (m_filetype == ctr::FieldFileType::EXODUSII) {
      // Create ExodusII writer
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
    }

  }
}

void
MeshWriter::writeMeta( uint64_t itr,
                       int chareid,
                       Centering centering,
                       const std::vector< std::string >& names )
// *****************************************************************************
//  Output metadata (field names) to file
//! \param[in] itr Iteration count since a new mesh (New mesh in this context
//!   means, either the mesh is moved and/or its topology has changed
//! \param[in] chareid The chare id the write-to-file request is coming from
//! \param[in] centering The centering that will be associated to the field data
//!   to be output when writeFields is called next
//! \param[in] names Names of fields to be output in next call to writeFields()
// *****************************************************************************
{
  if (!m_benchmark) {

    auto f = filename( itr, chareid );
  
    // Write field names
    #ifdef HAS_ROOT
    if (m_filetype == ctr::FieldFileType::ROOT) {
   
      // Create ROOT writer
      RootMeshWriter rmw( f, 1 );
      // Write node field names
      rmw.writeNodeVarNames( names );
  
    } else
    #endif
    if (m_filetype == ctr::FieldFileType::EXODUSII) {
  
      // Create ExodusII writer
      ExodusIIMeshWriter ew( f, ExoWriter::OPEN );
  
      // Write field names
      if (centering == Centering::ELEM)
        ew.writeElemVarNames( names );
      else if (centering == Centering::NODE)
        ew.writeNodeVarNames( names );
  
    }

  }
}

void
MeshWriter::writeFields( uint64_t itr,
                         uint64_t itf,
                         tk::real time,
                         int chareid,
                         Centering centering,
                         const std::vector< std::vector< tk::real > >& fields )
// *****************************************************************************
//  Output field data to file
//! \param[in] itr Iteration count since a new mesh (New mesh in this context
//!   means, either the mesh is moved and/or its topology has changed
//! \param[in] itf Field output iteration count
//! \param[in] time Physical time this at this field output dump
//! \param[in] chareid The chare id the write-to-file request is coming from
//! \param[in] centering The centering that will be associated to the field data
//!   to be output when writeFields is called next
//! \param[in] fields Field data to output to file
// *****************************************************************************
{
  if (!m_benchmark) {

std::cout << "wr: " << thisIndex << ", " << CkMyNode() << ", " << CkMyPe() << ": " << time << std::endl;

    auto f = filename( itr, chareid );
  
    // Write field names
    #ifdef HAS_ROOT
    if (m_filetype == ctr::FieldFileType::ROOT) {
   
      // Create ROOT writer
      RootMeshWriter rw( f, 1 );
      // Write time stamp
      rw.writeTimeStamp( itf, time );
      // Write node fields
      int varid = 0;
      for (const auto& v : fields) rw.writeNodeScalar( itf, ++varid, v );
  
    } else
    #endif
    if (m_filetype == ctr::FieldFileType::EXODUSII) {
  
      // Create ExodusII writer
      ExodusIIMeshWriter ew( f, ExoWriter::OPEN );
      // Write time stamp
      ew.writeTimeStamp( itf, time );
  
      // Write fields
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

std::string
MeshWriter::filename( uint64_t itr, int chareid ) const
// *****************************************************************************
//  Compute filename
//! \param[in] itr Iteration count since a new mesh (New mesh in this context
//!   means, either the mesh is moved and/or its topology has changed
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
  return m_outputBasefilename + ".e-s"
         + '.' + std::to_string( itr )          // create new file if new mesh
         + '.' + std::to_string( m_nchare )     // total number of workers
         + '.' + std::to_string( chareid )      // new file per worker
         #ifdef HAS_ROOT
         + (m_filetype == ctr::FieldFileType::ROOT ? ".root" : "")
         #endif
         ;
}

#include "NoWarning/meshwriter.def.h"
