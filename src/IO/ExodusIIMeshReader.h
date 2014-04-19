//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.h
  \author    J. Bakosi
  \date      Sat 19 Apr 2014 08:05:33 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader
*/
//******************************************************************************
#ifndef ExodusIIMeshReader_h
#define ExodusIIMeshReader_h

#include <string>

#include <UnsMesh.h>

namespace quinoa {

//! ExodusIIMeshReader
class ExodusIIMeshReader {

  public:
    //! Constructor
    explicit ExodusIIMeshReader( const std::string& filename,
                                 UnsMesh& mesh,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshReader();

    //! Read ExodusII mesh to file
    void read();

  private:
    //! Don't permit copy constructor
    ExodusIIMeshReader(const ExodusIIMeshReader&) = delete;
    //! Don't permit copy assigment
    ExodusIIMeshReader& operator=(const ExodusIIMeshReader&) = delete;
    //! Don't permit move constructor
    ExodusIIMeshReader(ExodusIIMeshReader&&) = delete;
    //! Don't permit move assigment
    ExodusIIMeshReader& operator=(ExodusIIMeshReader&&) = delete;

    //! Read header
    void readHeader();

    //! Read nodes
    void readNodes();

    //! Read elements
    void readElements();

    //! Read element block
    void readElemBlock();

    const std::string m_filename;          //!< File name

    UnsMesh& m_mesh;                       //!< Mesh object
    int m_inFile;                          //!< ExodusII file handle
    int m_neblk;                           //!< Number of element blocks
    int m_nnode;                           //!< Number of nodes
};

} // quinoa::

#endif // ExodusIIMeshReader_h
