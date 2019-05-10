// *****************************************************************************
/*!
  \file      src/IO/FileConvWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "NoWarning/TFile.hpp"
#include "NoWarning/TTree.hpp"
#include "NoWarning/TGraph2D.hpp"
#include "ExodusIIMeshWriter.hpp"
#include "Types.hpp"

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

    std::size_t m_nodal_size;

    std::unique_ptr< tk::ExodusIIMeshWriter > m_emw = nullptr;
    const std::string m_file_root;         //!< Root File name
    const std::string m_file_exodus;       //!< ExodusII File name

    //! Root file handle
    std::unique_ptr< TFile > m_infile = nullptr;
    //! ROOT tree handle
    TTree* m_tree_local = nullptr;

    //! Write the timestamp and the variables data
    void writeHeader();
    void writeCoordinates();
    void writeConnectivity();
    void writeVarNames();
    void writeData();

};

} // tk::

#endif // FileConvWriter_h
