// *****************************************************************************
/*!
  \file      src/IO/DiffWriterFiles.C
  \author    A. Pakki
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Comparator Writer Files
  \details   Compare the tolerance (absolute and relative) for files written by
     RootMeshWriter and ExodusIIMeshWriter classes. Right now, we are using the
     default absolute tolerance of 10e-6 and relative of 10e-9.
*/
// *****************************************************************************

#include <cmath>
#include <string>
#include <numeric>
#include <unordered_map>

#include "NoWarning/exodusII.h"
#include "NoWarning/TFile.h"
#include "NoWarning/TTree.h"

#include "DiffWriterFiles.h"
#include "ContainerUtil.h"
#include "Exception.h"

using tk::DiffWriterFiles;

DiffWriterFiles::DiffWriterFiles( const std::string& file_root,
                                  const std::string& file_exodus) :
  m_file_root( file_root ),
  m_file_exodus( file_exodus )
// *****************************************************************************
//  Constructor: open Exodus II file, Root file
//! \param[in] filename File to open as ExodusII file
//! \param[inout] mesh Unstructured mesh object to load the data to
//! \param[in] cpuwordsize Set CPU word size, see ExodusII documentation
//! \param[in] iowordsize Set I/O word size, see ExodusII documentation
//! \author A. Pakki
// *****************************************************************************
{

  emw = new tk::ExodusIIMeshWriter( file_exodus.c_str(), 
				    tk::ExoWriter::CREATE );

  m_infile = new TFile( file_root.c_str() );

  // If the ROOT file can be opened, fetch the TTree reference  
  if( m_infile == 0 ){
    std::cerr << "Unable to open the file " << file_root.c_str() << std::endl;
  } else {
    // the tree is named as ctree in the RootMeshWriter
    tree_local = (TTree*) m_infile->Get( "ctree" );
  }

}

DiffWriterFiles::~DiffWriterFiles() noexcept
// *****************************************************************************
//  Destructor
//! \author A. Pakki 
// *****************************************************************************
{
  // close the ROOT file
  m_infile->Close(); 

  delete emw;

}

void
DiffWriterFiles::readRootVar( uint64_t timestep, int varid,  
		  std::vector<double>** root_varvec )
//*****************************************************************************
//  Retrieve the variables from the ROOT writer
//  in[0]  - timestep
//  in[1]  - variable_id
//  out[0] - vector of values at the "timestep" time.
//  Author A. Pakki
//*****************************************************************************
{

  std::string branch_var = "branch_" + std::to_string(timestep) + "_field_"
			      + std::to_string(varid);

  tree_local->SetBranchAddress( branch_var.c_str(), root_varvec );
  tree_local->GetEntry( 0 );

}

void
DiffWriterFiles::convertFiles() 
//*****************************************************************************
//  Convert the input file [ROOT] layout to output file [ExodusII] layout
//  Author A. Pakki
//*****************************************************************************
{

  writeHeader();
  writeCoordinates();
  writeVarNames();
  //writeData();

}

void
DiffWriterFiles::writeHeader() 
//*****************************************************************************
//  Write the Header details
//  Author A. Pakki
//*****************************************************************************
{

  int coordinates, connectivity;
  /*
  tree_local->SetBranchAddress( "trian", &connectivity );
  tree_local->SetBranchAddress( "coord", &coordinates );
  
  tree_local->GetEntry( 0 );
  */
  emw->writeHeaderObject( "Data copied from ROOT",
                 3, 
		 //coordinates,
		 //(connectivity / 4),
		 1201, 3643, 
		 1, 0, 0 );
}

void
DiffWriterFiles::writeCoordinates()
//*****************************************************************************
//  Write the Coordinates (x,y,z) data
//  Author A. Pakki
//*****************************************************************************
{

  std::vector< tk::real > *mx = nullptr;
  std::vector< tk::real > *my = nullptr;
  std::vector< tk::real > *mz = nullptr;

  tree_local->SetBranchAddress( "x_coord", &mx );
  tree_local->SetBranchAddress( "y_coord", &my );
  tree_local->SetBranchAddress( "z_coord", &mz );
  tree_local->GetEntry( 0 );

  // write to ExodusII
  emw->writeNodesObject( *mx, *my, *mz );
}

void
DiffWriterFiles::writeVarNames()
//*****************************************************************************
//  Write the Variables names
//  Author A. Pakki
//*****************************************************************************
{
  std::vector<std::string> *var_copy = nullptr;
  
  tree_local->SetBranchAddress( "variables", &var_copy );
  tree_local->GetEntry(0);

  emw->writeNodeVarNames( *var_copy );

  variables = (*var_copy).size();
}

void
DiffWriterFiles::writeData()
//*****************************************************************************
//  Write the Timestamp and Variables from Input to Output
//  Author A. Pakki
//*****************************************************************************
{

  int timestep = 1; 
  while( true ) {
    double dt;

    std::string time_branch =  "time_branch_" + std::to_string(timestep);
    if( tree_local->GetBranch( time_branch.c_str() ) == nullptr )
      break;

    for (int var_id = 1; var_id <= variables; var_id++ ) {
      std::string branch_var = "branch_" + std::to_string(timestep) + "_field_"
			      + std::to_string(var_id);

      if( tree_local->GetBranch( branch_var.c_str() ) == nullptr )
	break;
      else {

	std::vector< double> *var_fields   = nullptr;

	tree_local->SetBranchAddress( branch_var.c_str(), &var_fields );
	tree_local->SetBranchAddress(time_branch.c_str(), &dt );
	tree_local->GetEntry(0);
	
	emw->writeNodeScalar( timestep, var_id, *var_fields );
	//tree_local->ResetBranchAddresses();
      }	

    } // End the variables loop.
    
    // Write the timestamp variable, once per timestep
    emw->writeTimeStamp( timestep, dt );
    timestep++;

  } // break after all the variables for all the timesteps are written
}
