//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshReader.h
  \author    J. Bakosi
  \date      Wed Apr 23 08:42:19 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     ExodusII mesh reader
  \details   ExodusII mesh reader
*/
//******************************************************************************
#ifndef ExodusIIMeshReader_h
#define ExodusIIMeshReader_h

#include <string>

#include <UnsMesh.h>
#include <Reader.h>

namespace quinoa {

//! ExodusIIMeshReader
class ExodusIIMeshReader : public tk::Reader {

  public:
    //! Constructor
    explicit ExodusIIMeshReader( const std::string& filename,
                                 UnsMesh& mesh,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshReader() noexcept override;

    //! Read ExodusII mesh to file
    void read() override;

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
