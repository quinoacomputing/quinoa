// *****************************************************************************
/*!
  \file      src/IO/FileConvWriter.C
  \author    A. Pakki
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Converter Writer Files
  \details   Convert the input file written by RootMeshWriter into ExodusII 
             layout. We make use of the ExodusIIMeshWriter class to write the 
	     output file.
*/
// *****************************************************************************

#include <cmath>
#include <string>
#include <numeric>

#include "Make_unique.h"

#include "NoWarning/exodusII.h"
#include "NoWarning/TFile.h"
#include "NoWarning/TTree.h"

#include "FileConvWriter.h"
#include "ContainerUtil.h"
#include "Exception.h"

using tk::FileConvWriter;

FileConvWriter::FileConvWriter( const std::string& file_root,
                                const std::string& file_exodus) :
  m_file_root( file_root ),
  m_file_exodus( file_exodus )
// *****************************************************************************
//  Constructor: open Exodus II file, Root file
//! \param[in] file_root File to open as ROOT file
//! \param[in] file_exodus File to be opened as ExodusII file to be written
//! \author A. Pakki
// *****************************************************************************
{

  m_emw = tk::make_unique< tk::ExodusIIMeshWriter >
                         ( file_exodus.c_str(), tk::ExoWriter::CREATE );

  m_infile = tk::make_unique< TFile >( file_root.c_str() );

  // If the ROOT file can be opened, fetch the TTree reference  
  if( m_infile == nullptr ){
    std::cerr << "Unable to open the file " << file_root.c_str() << std::endl;
  } else {
    // the tree is named as ctree in the RootMeshWriter
    m_tree_local = static_cast< TTree* >( m_infile->Get( "ctree" ) );
  }

}

FileConvWriter::~FileConvWriter() noexcept
// *****************************************************************************
//  Destructor
//! \author A. Pakki 
// *****************************************************************************
{
  // close the ROOT file
  m_infile->Close(); 
}

void
FileConvWriter::convertFiles() 
//*****************************************************************************
//  Convert the input file [ROOT] layout to output file [ExodusII] layout
//  Author A. Pakki
//*****************************************************************************
{

  writeHeader();
  writeVarNames();
  writeConnectivity();
  writeCoordinates();
  writeData();

}

void
FileConvWriter::writeHeader() 
//*****************************************************************************
//  Write the Header details
//  Author A. Pakki
//*****************************************************************************
{

  int coord, connect;

  m_tree_local->SetBranchAddress( "trian", &connect );
  m_tree_local->SetBranchAddress( "coord", &coord );

  m_tree_local->GetEntry( 0 );
  m_emw->writeHeader( "Data copied from ROOT",
                 3, 
		 coord,
		 (connect / 4),
		 1, 0, 0 );
  
  m_tree_local->ResetBranchAddresses();

}

void
FileConvWriter::writeCoordinates()
//*****************************************************************************
//  Write the Coordinates (x,y,z) data
//  Author A. Pakki
//*****************************************************************************
{

  std::vector< tk::real > *mx = nullptr;
  std::vector< tk::real > *my = nullptr;
  std::vector< tk::real > *mz = nullptr;

  m_tree_local->SetBranchAddress( "x_coord", &mx );
  m_tree_local->SetBranchAddress( "y_coord", &my );
  m_tree_local->SetBranchAddress( "z_coord", &mz );

  m_tree_local->GetEntry( 0 );

  // write to ExodusII
  m_emw->writeNodes( *mx, *my, *mz );
  m_tree_local->ResetBranchAddresses();

}

void
FileConvWriter::writeConnectivity()
//*****************************************************************************
//  Write the Connectivity between the coordinates.
//  Author A. Pakki
//*****************************************************************************
{

  // param[0] - elclass, refer ExodusIIMeshWriter.C
  // param[1] - vertices, 4 is the number of vertices of Tetrahedron
  // param[2] - string literal TETRAHEDRA/TRIANGLES
  int elclass = 0;
  std::vector<std::size_t> *tets_number = nullptr;

  m_tree_local->SetBranchAddress( "tetconnect", &tets_number );
  m_tree_local->GetEntry( 0 );

  m_emw->writeElemBlock( elclass, 4, "TETRAHEDRA", *tets_number );
  m_tree_local->ResetBranchAddresses();

}

void
FileConvWriter::writeVarNames()
//*****************************************************************************
//  Write the Variables names
//  Author A. Pakki
//*****************************************************************************
{
  std::vector<std::string> *var_copy = nullptr;
  
  m_tree_local->SetBranchAddress( "variables", &var_copy );
  m_tree_local->GetEntry(0);

  m_emw->writeNodeVarNames( *var_copy );
  m_nodal_size = var_copy->size();
  m_tree_local->ResetBranchAddresses();

}

void
FileConvWriter::writeData()
//*****************************************************************************
//  Write the Timestamp and Variables from Input to Output
//  Author A. Pakki
//*****************************************************************************
{

  std::size_t timestep = 1; 
  while( true ) {
    double dt = 0;

    std::string time_branch =  "time_branch_" + std::to_string(timestep);
    if( m_tree_local->GetBranch( time_branch.c_str() ) == nullptr )
      break;

    for (std::size_t var_id=1; var_id<=m_nodal_size; ++var_id) {
      std::string branch_var = "branch_" + std::to_string(timestep) + "_field_"
			      + std::to_string(var_id);

      if( m_tree_local->GetBranch( branch_var.c_str() ) == nullptr )
	break;
      else {

	std::vector< double> *var_fields   = nullptr;

	m_tree_local->SetBranchAddress( branch_var.c_str(), &var_fields );
	m_tree_local->SetBranchAddress(time_branch.c_str(), &dt );
	m_tree_local->GetEntry(0);
	
	m_emw->writeNodeScalar( timestep,
                                static_cast<int>(var_id),
                                *var_fields );
      }	

    } // End the variables loop.
    
    // Write the timestamp variable, once per timestep
    m_emw->writeTimeStamp( timestep, dt );
    ++timestep;

  } // break after all the variables for all the timesteps are written
  m_tree_local->ResetBranchAddresses();

}

