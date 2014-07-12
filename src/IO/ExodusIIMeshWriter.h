//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.h
  \author    J. Bakosi
  \date      Wed Apr 23 11:12:23 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     ExodusII mesh writer
  \details   ExodusII mesh writer
*/
//******************************************************************************
#ifndef ExodusIIMeshWriter_h
#define ExodusIIMeshWriter_h

#include <string>

#include <UnsMesh.h>
#include <Writer.h>

namespace quinoa {

//! ExodusIIMeshWriter
class ExodusIIMeshWriter : public tk::Writer {

  public:
    //! Constructor
    explicit ExodusIIMeshWriter( const std::string& filename,
                                 UnsMesh& mesh,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshWriter() noexcept override;

    //! Write ExodusII mesh to file
    void write() override;

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
