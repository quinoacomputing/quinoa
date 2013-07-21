//******************************************************************************
/*!
  \file      src/IO/SiloWriter.h
  \author    J. Bakosi
  \date      Sat 20 Jul 2013 06:39:02 PM MDT
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

namespace Quinoa {

class STLMesh;

//! SiloWriter
class SiloWriter {

  public:
    //! Constructor: Acquire glob file handle
    explicit SiloWriter(const std::string& filename,
                        STLMesh* const mesh,
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

    const std::string m_filename;     //!< Silo filename
    STLMesh* const m_mesh;            //!< Mesh object pointer

    DBfile* m_dbfile;                 //!< Silo DB file
};

} // namespace Quinoa

#endif // SiloWriter_h
