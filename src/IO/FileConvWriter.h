// *****************************************************************************
/*!
  \file      src/IO/FileConvWriter.h
  \author    A. Pakki 
  \copyright 2012-2015, Aditya Pakki, 2016, Los Alamos National Security, LLC.
  \brief     Convert the input files to output file format  
  \details   Convert the input file in ROOT format to the ExodusII layout
	     The next step post the output would be to utilize exodiff utility.
*/
// *****************************************************************************
#ifndef FileConvWriter_h
#define FileConvWriter_h

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "NoWarning/TFile.h"
#include "NoWarning/TTree.h"
#include "NoWarning/TGraph2D.h"
#include "ExodusIIMeshWriter.h"

#include "Types.h"

namespace tk {


class FileConvWriter {

  public:
    //! Constructor: Open ROOT, ExodusII files
    explicit FileConvWriter( const std::string& file_root, 
			      const std::string& file_exodus );

    //! Destructor
    ~FileConvWriter() noexcept;
  
    //! Convert the files
    void convertFiles();

  private:

    int nodal_size;

    tk::ExodusIIMeshWriter *emw = nullptr;
    const std::string m_file_root;         // Root File name
    const std::string m_file_exodus;       // ExodusII File name

    TFile *m_infile = nullptr;                // Root File handle
    TTree *tree_local = nullptr;           // Get the tree handle from ROOT

    //! Write the timestamp and the variables data
    void writeHeader();
    void writeCoordinates();
    void writeConnectivity();
    void writeVarNames();
    void writeData();

};

} // tk::

#endif // FileConvWriter_h
