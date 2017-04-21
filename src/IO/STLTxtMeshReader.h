// *****************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     ASCII STL (STereoLithography) reader class declaration
  \details   ASCII STL (STereoLithographu) reader class declaration.
*/
// *****************************************************************************
#ifndef STLTxtMeshReader_h
#define STLTxtMeshReader_h

#include <iostream>

#include "Reader.h"
#include "Exception.h"

class STLMesh;

namespace tk {

//! \brief STLTxtMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a file in
//!   ASCII STL format.
class STLTxtMeshReader : public Reader {

  public:
    //! Constructor
    explicit STLTxtMeshReader( const std::string filename, STLMesh& mesh );

    //! Read ASCII STL mesh
    void readMesh();

  private:
    //! \brief ASCII STL keyword with operator>> redefined to do error checking
    //!    without contaminating client-code
    //! \author J. Bakosi
    struct STLKeyword {
      std::string read;                 //!< Keyword read in from input
      const std::string correct;        //!< Keyword that should be read in

      //! Initializer constructor
      explicit STLKeyword( const std::string& corr ) noexcept :
        read(), correct(corr) {}

      //! Operator >> for reading a keyword and hande error
      friend std::ifstream& operator>> (std::ifstream& is, STLKeyword& kw) {
        is >> kw.read;
        ErrChk( kw.read == kw.correct,
                "Corruption in ASCII STL file while parsing keyword '" +
                kw.read + "', should be '" + kw.correct + "'" );
        return is;
      }
    };

    //! Read (or count vertices in) ASCII STL mesh
    std::size_t readFacets( const bool store,
                            tk::real* const x = nullptr,
                            tk::real* const y = nullptr,
                            tk::real* const z = nullptr );

    const bool STORE = true;                 //!< Indicator to store facets
    const bool COUNT = false;                //!< Indicator to only count facets
    STLMesh& m_mesh;                         //!< Mesh
};

} // tk::

#endif // STLTxtMeshReader_h
