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

#include "Types.h"

namespace tk {

class UnsMesh;

class DiffWriterFiles {

  public:
    //! Constructor: Open ROOT, ExodusII files
    explicit DiffWriterFiles( const std::string& file_root, 
			      const std::string& file_exodus,
                              int cpuwordsize = sizeof(double),
                              int iowordsize = sizeof(double) );

    //! Destructor
    ~DiffWriterFiles() noexcept;
  
    void computeDifferences();

  private:

    const std::string m_file_root;         // Root File name
    const std::string m_file_exodus;       // ExodusII File name

    int m_inFile;                       //!< ExodusII file handle
    TFile *rfile = nullptr;             // Root File handle
    TTree *tree_local = nullptr;      // Get the tree handle in the file

    int variables;			// stored variables
    uint64_t timesteps;			// stored timesteps
    int nodal_size;                     // number of nodes

    double abs_tolerance = 10e-6;
    double rel_tolerance = 10e-9;

    // read ExodusII file timesteps
    void readExodusIIVar( uint64_t timestep, int varid,
			  tk::real* exodus_varvec );

    // read ROOT file timesteps
    void readRootVar( uint64_t timestep, int varid, 
		      std::vector<double>** root_varvec );

};

} // tk::

#endif // DiffWriterFiles_h
