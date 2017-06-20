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

#include <cstdint>
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
#include "UnsMesh.h"
#include "Reorder.h"

using tk::DiffWriterFiles;

DiffWriterFiles::DiffWriterFiles( const std::string& file_root,
                                  const std::string& file_exodus,
				  int cpuwordsize,
                                  int iowordsize ) :
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
  float version;

  m_inFile = ex_open( file_exodus.c_str(), EX_READ, &cpuwordsize, &iowordsize,
                      &version );

  ErrChk( m_inFile > 0, "Failed to open ExodusII file: " + file_exodus );

  rfile = new TFile( file_root.c_str() );

  // If the ROOT file can be opened, fetch the TTree reference  
  if( rfile == 0 ){
    std::cerr << "Unable to open the file " << file_root.c_str() << std::endl;
  } else {
    // the tree is named as ctree in the RootMeshWriter
    tree_local = (TTree*) rfile->Get( "ctree" );
  }

}

DiffWriterFiles::~DiffWriterFiles() noexcept
// *****************************************************************************
//  Destructor
//! \author A. Pakki 
// *****************************************************************************
{
  if ( ex_close(m_inFile) < 0 )
    printf( ">>> WARNING: Failed to close ExodusII file: %s\n",
            m_file_exodus.c_str() );

  // close the ROOT file
  rfile->Close(); 
}


void
DiffWriterFiles::computeDifferences() 
// *****************************************************************************
//  Perform the comparison of the values for the variables
//! \author A. Pakki
// *****************************************************************************
{
  // Fetch the number of variables stored in the File.
  ErrChk(
    ex_get_variable_param( m_inFile, EX_NODE_BLOCK,
                           &variables )  == 0,
    "Failed to write nodal output variable parameters to ExodusII file: " +
    m_file_exodus );

  // fetch the variables from ROOT layout.
  std::vector<std::string> *var_copy = nullptr;
  
  tree_local->SetBranchAddress( "variables", &var_copy );

  tree_local->GetEntry(0);
  Assert( variables == var_copy->size(), "The variables count are not same" );

  //Fetch the number of timesteps
  timesteps =  ex_inquire_int( m_inFile, EX_INQ_TIME );
  Assert( timesteps > 0, "No timesteps found to query in " + m_file_exodus );

  //Fetch the nodal points
  nodal_size = ex_inquire_int( m_inFile, EX_INQ_NODES );
  Assert( nodal_size > 0, "No nodes found to query in " + m_file_exodus );

  for( uint64_t timeItr = 1; timeItr <= timesteps; timeItr++ ) {

    for (int varItr = 1; varItr <= variables; varItr++ ) {

      tk::real* exodus_varvec = new tk::real[nodal_size];
      std::vector<double>* root_varvec = nullptr;
 
      // ExodusII is a 1 based indexing and ROOT is 0 based indexing
      readExodusIIVar( timeItr, varItr, exodus_varvec );
      
      readRootVar( timeItr, varItr, &root_varvec );

      Assert( nodal_size == root_varvec->size(), "Nodal points are not same" );

      double difference = 0.0;
      for( int i = 0 ; i < nodal_size; i++) 
	difference += abs( static_cast<double>( exodus_varvec[i] ) - 
		    (*root_varvec)[i] );
      
      Assert( difference < abs_tolerance, "Difference exceeds tolerance" );
      delete exodus_varvec;
    }

  }
  
}

void
DiffWriterFiles::readExodusIIVar( uint64_t timestep, int varid, 
		  tk::real* exodus_varvec ) 
//*****************************************************************************
//  Retrieve the variables from the ExodusII writer
//  in[0]  - timestep
//  in[1]  - variable_id
//  out[0] - vector of values at the "timestep" time.
//  Author A. Pakki
//*****************************************************************************
{
  
  ErrChk( ex_get_var( m_inFile,
                      static_cast< int >( timestep ),
                      EX_NODE_BLOCK,
                      varid,
                      1,
                      nodal_size,
                      exodus_varvec ) == 0,
          "Failed to read node scalar from ExodusII file: " + m_file_exodus );

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

