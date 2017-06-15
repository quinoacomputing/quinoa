// *****************************************************************************
/*!
  \file      src/IO/RootMeshWriter.C
  \author    A. Pakki
  \copyright 2012-2015, Jozsef Pakki, 2016, Los Alamos National Security, LLC.
  \brief     Root mesh-based data writer
  \details   Root mesh-based data writer class definition.
*/
// *****************************************************************************

#include <algorithm>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>

#include "RootMeshWriter.h"
#include "Exception.h"
#include "UnsMesh.h"

using tk::RootMeshWriter;

RootMeshWriter::RootMeshWriter( const std::string filename, int option ) :
  m_filename( filename )
// *****************************************************************************
//  Constructor: create/open Root file
//! \param[in] "filename" File to open as Root file
//! \param[in] mode Root writer constructor mode: ExoWriter::CREATE for
//!   creating a new file, ExoWriter::OPEN for opening an existing file for
//!   appending
//! \author A. Pakki
// *****************************************************************************
{
  if (option == 0 ) {

    rfile = new TFile(filename.c_str(), "RECREATE" );
    tree_connect = new TTree ( "ctree", "store the connectivity" );

  } else if (option == 1) {

      rfile = TFile::Open(filename.c_str(), "UPDATE" );
      tree_connect = (TTree*) rfile->Get( "ctree" );

  } else Throw( "Root Mesh modes not supported" );
}

RootMeshWriter::~RootMeshWriter() noexcept
// *****************************************************************************
//  Destructor
//! \author A. Pakki
// *****************************************************************************
{
  if (rfile) {
    rfile->Write();
    rfile->Close();
  }
}

void
RootMeshWriter::writeMesh( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write Root mesh file
//! \param[in] mesh Unstructured mesh object
//! \author A. Pakki
// *****************************************************************************
{
  writeHeader( mesh );
  writeNodes( mesh );
  writeElements( mesh );
}

void
RootMeshWriter::writeHeader( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write Root header
//! \param[in] mesh Unstructured mesh object
//! \author A. Pakki
// *****************************************************************************
{

  csobject = new mesh_data(mesh.size(), ( mesh.tetinpoel().size() + 
		  mesh.triinpoel().size() ) );

}

void
RootMeshWriter::writeNodes( const UnsMesh& mesh )  const
// *****************************************************************************
//  Write node coordinates to Root file
//! \param[in] mesh Unstructured mesh object
//! \author A. Pakki
// *****************************************************************************
{

  // the file requires the vertices and the number of triangles
  // 4 triangles per tetrahedron and mesh.tetinpoel() stores 4 
  // vertices per tet in the vector (# vertices = # of triangles)

  tree_connect->Branch( "coord", &csobject->coordinates, "coordinates/I" );
  tree_connect->Branch( "trian", &csobject->triangles, "triangles/I" );
  
  tree_connect->Branch( "x_coord", &csobject->mx_root );
  tree_connect->Branch( "y_coord", &csobject->my_root );
  tree_connect->Branch( "z_coord", &csobject->mz_root );

  for ( int i = 0 ; i < csobject->coordinates; i++ ) {
    csobject->mx_root.push_back( mesh.x()[i] );
    csobject->my_root.push_back( mesh.y()[i] );
    csobject->mz_root.push_back( mesh.z()[i] );
  }
}

void
RootMeshWriter::writeElements( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write element connectivity to Root file
//! \param[in] mesh Unstructured mesh object
//! \author A. Pakki
// *****************************************************************************
{
  int elclass = 0;

  writeElemBlock( elclass, mesh.triinpoel() );
  writeElemBlock( elclass, mesh.tetinpoel() );
}

void
RootMeshWriter::writeElemBlock( int& elclass,
                                const std::vector< std::size_t >& inpoel )
const
// *****************************************************************************
//  Write element block to ROOT file
//! \param[inout] elclass Count element class ids in file
//! \param[in] inpoel Element connectivity.
//! \author A. Pakki
// *****************************************************************************
{
  if (inpoel.empty()) return;

  // increase number of element classes in file
  ++elclass;

  // Make sure element connectivity starts with zero
  Assert( *std::minmax_element( begin(inpoel), end(inpoel) ).first == 0,
          "node ids should start from zero" );
  
  // create a branch for storing the tetrahedrons  
  tree_connect->Branch( "connect", &csobject->connectivity );
  for ( auto itr = inpoel.begin(); itr != inpoel.end(); itr++ )
    csobject->connectivity.push_back( *itr );

  tree_connect->Fill();
  tree_connect->Write(); 
}

void
RootMeshWriter::writeNodeVarNames( const std::vector< std::string >& nv )
const
// *****************************************************************************
//  Write the names of nodal output variables to ROOT file
//! \param[in] nv Nodal variable names
//! \author A. Pakki
// *****************************************************************************
{
  std::vector < std::string > nv_copy;
  TBranch *branch = tree_connect->Branch( "variables", &nv_copy );

  for ( auto itr = nv.begin(); itr != nv.end(); itr++)
    nv_copy.push_back( *itr );

  branch->Fill();
  tree_connect->Write();

}

void
RootMeshWriter::writeTimeStamp( uint64_t it, tk::real time ) const
// *****************************************************************************
//  Write time stamp to ROOT file
//! \param[in] it Iteration number
//! \param[in] time Time
//! \author A. Pakki
// *****************************************************************************
{
  
  // create a friend tree to the main tree
  std::string tree = "tf_var_" + std::to_string(it);
  friendTree = new TTree( tree.c_str(), "friend trees" );

  // create a branch for the time step
  tk::real elapsed_time = time;
  std::string branch_time = "time_branch_" + std::to_string(it);
  friendTree->Branch( branch_time.c_str(), &elapsed_time, "elapsed_time/D" );

  friendTree->Fill();
  tree_connect->AddFriend( tree.c_str() );

}

void
RootMeshWriter::writeNodeScalar( uint64_t it,
                                     int varid,
                                     const std::vector< tk::real >& var ) const
// *****************************************************************************
//  Write node scalar field to ROOT file
//! \param[in] it Iteration number
//! \param[in] varid Variable id
//! \param[in] var Vector of variable to output
//! \author A. Pakki
// *****************************************************************************
{

  std::vector< tk::real > vec_copy;

  std::string branch_field = "branch_" + std::to_string(it) + "_field_" + 
			      std::to_string(varid);
  //create a Branch
  friendTree->Branch( branch_field.c_str(), &vec_copy );
  
  //copy the values
  for( auto itr = var.begin(); itr != var.end(); itr++ )
    vec_copy.push_back( *itr );

  friendTree->Fill();
  friendTree->Write();

}

