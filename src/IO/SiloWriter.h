//******************************************************************************
/*!
  \file      src/IO/SiloWriter.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 10:07:15 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Silo writer declaration.
  \details   Silo writer declaration. This class currently only supports writing
    an STL triangulation into a Silo file. See also
    https://wci.llnl.gov/codes/silo.
*/
//******************************************************************************
#ifndef SiloWriter_h
#define SiloWriter_h

#include <string>

#include <silo.h>

namespace quinoa {

//! Silo error handler function type
typedef void (*SiloErrorHandler)( char* );

//! Silo error handler
void SiloError( char* msg );

class STLMesh;

//! \brief SiloWriter
//! \details Mesh reader class facilitating reading a mesh from a file in
//!   Silo format. See also https://wci.llnl.gov/codes/silo.
class SiloWriter {

  public:
    //! Constructor
    explicit SiloWriter( const std::string& filename,
                         const STLMesh& mesh,
                         int errLevel );
 
    //! Destructor: Release glob file handle
    ~SiloWriter() noexcept;

    //! Write out silo file
    void write();

  private:
    const std::string m_filename;       //!< Silo filename
    const STLMesh& m_mesh;              //!< Mesh

    SiloErrorHandler m_errFunc;         //!< Silo error handler function ptr
    int m_errLevel;                     //!< Silo error reporting level
    DBfile* m_dbfile;                   //!< Silo DB file
};

} // quinoa::

#endif // SiloWriter_h
