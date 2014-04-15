//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.h
  \author    J. Bakosi
  \date      Tue Apr 15 07:52:14 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ExodusII mesh writer
  \details   ExodusII mesh writer
*/
//******************************************************************************
#ifndef ExodusIIMeshWriter_h
#define ExodusIIMeshWriter_h

#include <string>

#include <UnsMesh.h>

namespace quinoa {

//! ExodusIIMeshWriter
class ExodusIIMeshWriter {

  public:
    //! Constructor
    explicit ExodusIIMeshWriter( const std::string& filename,
                                 UnsMesh& mesh,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshWriter();

    //! Write ExodusII mesh to file
    void write();

  private:
    //! Don't permit copy constructor
    ExodusIIMeshWriter(const ExodusIIMeshWriter&) = delete;
    //! Don't permit copy assigment
    ExodusIIMeshWriter& operator=(const ExodusIIMeshWriter&) = delete;
    //! Don't permit move constructor
    ExodusIIMeshWriter(ExodusIIMeshWriter&&) = delete;
    //! Don't permit move assigment
    ExodusIIMeshWriter& operator=(ExodusIIMeshWriter&&) = delete;

    //! Write header
    void writeHeader();

    //! Write nodes
    void writeNodes();

    //! Write elements
    void writeElements();

    //! Write element block
    void writeElemBlock( int elclass, int nnpe, const std::string& eltype,
                         std::vector< std::vector< int > >& inpoel );

    const std::string m_filename;          //!< File name

    UnsMesh& m_mesh;                       //!< Mesh object
    int m_outFile;                         //!< ExodusII file handle
};

} // quinoa::

#endif // ExodusIIMeshWriter_h
