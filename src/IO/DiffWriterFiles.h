// *****************************************************************************
/*!
  \file      src/IO/DiffWriterFiles.h
  \author    A. Pakki 
  \copyright 2012-2015, Aditya Pakki, 2016, Los Alamos National Security, LLC.
  \brief     Compare the output files written by two Writers
  \details   Compute the relative and absolute differences for variables
	     written using RootMeshWriter and ExodusII Mesh Writer.
*/
// *****************************************************************************
#ifndef DiffWriterFiles_h
#define DiffWriterFiles_h

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "NoWarning/TFile.h"
#include "NoWarning/TTree.h"
#include "NoWarning/TGraph2D.h"
#include "ExodusIIMeshWriter.h"

#include "Types.h"

namespace tk {


class DiffWriterFiles {

  public:
    //! Constructor: Open ROOT, ExodusII files
    explicit DiffWriterFiles( const std::string& file_root, 
			      const std::string& file_exodus );

    //! Destructor
    ~DiffWriterFiles() noexcept;
  
    //! Convert the files
    void convertFiles();

  private:

    int variables;

    tk::ExodusIIMeshWriter *emw = nullptr;
    const std::string m_file_root;         // Root File name
    const std::string m_file_exodus;       // ExodusII File name

    TFile *m_infile = nullptr;                // Root File handle
    TTree *tree_local = nullptr;           // Get the tree handle from ROOT

    // read ROOT file timesteps
    void readRootVar( uint64_t timestep, int varid, 
		      std::vector<double>** root_varvec );
    
    //! Write the timestamp and the variables data
    void writeHeader();
    void writeCoordinates();
    void writeVarNames();
    void writeData();

};

} // tk::

#endif // DiffWriterFiles_h
