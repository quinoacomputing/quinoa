//******************************************************************************
/*!
  \file      src/IO/SiloWriter.h
  \author    J. Bakosi
  \date      Fri Sep 27 15:07:28 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Silo (https://wci.llnl.gov/codes/silo) writer
  \details   Silo (https://wci.llnl.gov/codes/silo) writer
*/
//******************************************************************************
#ifndef SiloWriter_h
#define SiloWriter_h

#include <string>

#include <silo.h>

#include <QuinoaTypes.h>

namespace quinoa {

//! Silo error handler function type
typedef void (*SiloErrorHandler)(char*);

//! Silo error handler
void SiloError(char* msg);

class STLMesh;

//! SiloWriter
class SiloWriter {

  public:
    //! Constructor: Acquire glob file handle
    explicit SiloWriter(const std::string& filename,
                        STLMesh& mesh,
                        const int errLevel);
 
    //! Destructor: Release glob file handle
    ~SiloWriter() noexcept;

    //! Write out silo file
    void write();

  private:
    //! Don't permit copy constructor
    SiloWriter(const SiloWriter&) = delete;
    //! Don't permit copy assigment
    SiloWriter& operator=(const SiloWriter&) = delete;
    //! Don't permit move constructor
    SiloWriter(SiloWriter&&) = delete;
    //! Don't permit move assigment
    SiloWriter& operator=(SiloWriter&&) = delete;

    const std::string m_filename;       //!< Silo filename
    STLMesh& m_mesh;                    //!< Mesh

    SiloErrorHandler m_errFunc;         //!< Silo error handler function ptr
    int m_errLevel;                     //!< Silo error reporting level
    DBfile* m_dbfile;                   //!< Silo DB file
};

} // namespace quinoa

#endif // SiloWriter_h
