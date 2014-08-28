//******************************************************************************
/*!
  \file      src/Base/Writer.h
  \author    J. Bakosi
  \date      Thu 28 Aug 2014 11:51:40 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Writer base class declaration
  \details   Writer base class declaration
*/
//******************************************************************************
#ifndef Writer_h
#define Writer_h

#include <string>
#include <fstream>

#include <Exception.h>

namespace tk {

//! Writer base class
class Writer {

  public:
    //! Constructor: Acquire file handle
    explicit Writer(const std::string& filename);

    //! Destructor: Release file handle
    virtual ~Writer() noexcept;

    //! Write interface
    virtual void write() { Throw( "Writer::write() is a no-op" ); }

  protected:
    const std::string m_filename;          //!< File name
    mutable std::ofstream m_outFile;       //!< File output stream
};

} // tk::

#endif // Writer_h
