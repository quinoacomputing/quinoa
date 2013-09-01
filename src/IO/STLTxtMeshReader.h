//******************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.h
  \author    J. Bakosi
  \date      Sun 01 Sep 2013 02:22:01 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) reader class declaration
  \details   ASCII STL (STereoLithographu) reader class declaration
*/
//******************************************************************************
#ifndef STLTxtMeshReader_h
#define STLTxtMeshReader_h

#include <vector>
#include <iostream>

#include <Reader.h>

class STLMesh;

namespace quinoa {

//! STLTxtMeshReader : Reader
class STLTxtMeshReader : public Reader {

  public:
    //! Constructor
    explicit STLTxtMeshReader(const std::string filename,
                              STLMesh* const mesh) :
      Reader(filename),
      m_mesh(mesh) {
      Assert(m_mesh != nullptr, ExceptType::FATAL,
            "Uninitialized mesh object passed to STLTxtMeshReader constructor");
      // Set mesh name as filename modulo extension
      mesh->setName(filename.substr(0, filename.find_last_of(".")));
    }

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
      explicit STLKeyword(const std::string corr) noexcept : correct(corr) {}

      //! Operator >> for reading a keyword and error handling
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
    STLMesh* const m_mesh;                   //!< Mesh object pointer
};

} // namespace quinoa

#endif // STLTxtMeshReader_h
