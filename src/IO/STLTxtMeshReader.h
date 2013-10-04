//******************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.h
  \author    J. Bakosi
  \date      Thu 03 Oct 2013 09:00:44 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) reader class declaration
  \details   ASCII STL (STereoLithographu) reader class declaration
*/
//******************************************************************************
#ifndef STLTxtMeshReader_h
#define STLTxtMeshReader_h

#include <iostream>

#include <Reader.h>
#include <Exception.h>

class STLMesh;

namespace quinoa {

//! STLTxtMeshReader : Reader
class STLTxtMeshReader : public Reader {

  public:
    //! Constructor
    explicit STLTxtMeshReader(const std::string filename, STLMesh& mesh);

    //! Destructor, default compiler generated
    ~STLTxtMeshReader() noexcept override = default;

    //! Read ASCII STL mesh
    void read() override;

  private:
    //! Don't permit copy constructor
    STLTxtMeshReader(const STLTxtMeshReader&) = delete;
    //! Don't permit copy assigment
    STLTxtMeshReader& operator=(const STLTxtMeshReader&) = delete;
    //! Don't permit move constructor
    STLTxtMeshReader(STLTxtMeshReader&&) = delete;
    //! Don't permit move assigment
    STLTxtMeshReader& operator=(STLTxtMeshReader&&) = delete;

    //! ASCII STL keyword with operator>> redefined to do error checking
    struct STLKeyword {
      std::string read;                 //!< Keyword read in from input
      const std::string correct;        //!< Keyword that should be read in

      //! Initializer constructor
      explicit STLKeyword(const std::string& corr) noexcept : correct(corr) {}

      //! Operator >> for reading a keyword and hande error
      friend std::ifstream& operator>> (std::ifstream& is, STLKeyword& kw) {
        is >> kw.read;
        ErrChk(kw.read == kw.correct, ExceptType::FATAL,
               "Corruption in ASCII STL file while parsing keyword '" +
               kw.read + "', should be '" + kw.correct + "'");
        return is;
      }
    };

    //! Read (or count vertices in) ASCII STL mesh
    size_t readFacets(const bool store,
                      real* const x = nullptr,
                      real* const y = nullptr,
                      real* const z = nullptr);

    const bool STORE = true;                 //!< Indicator to store facets
    const bool COUNT = false;                //!< Indicator to only count facets
    STLMesh& m_mesh;                         //!< Mesh
};

} // namespace quinoa

#endif // STLTxtMeshReader_h
