//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.h
  \author    J. Bakosi
  \date      Sat 13 Jul 2013 08:03:18 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     STL (STereoLithography) class declaration
  \details   STL (STereoLithography) class declaration
*/
//******************************************************************************
#ifndef STLMesh_h
#define STLMesh_h

#include <QuinoaTypes.h>

namespace Quinoa {

//! STLMesh
class STLMesh {

  public:
    //! Constructor, default compiler generated
    explicit STLMesh() = default;

    //! Destructor, default compiler generated
    ~STLMesh() noexcept = default;

  private:
    //! Don't permit copy constructor
    STLMesh(const STLMesh&) = delete;
    //! Don't permit assigment constructor
    STLMesh& operator=(const STLMesh&) = delete;
    //! Don't permit move constructor
    STLMesh(STLMesh&&) = delete;
    //! Don't permit move assignment
    STLMesh& operator=(STLMesh&&) = delete;
};

} // namespace Quinoa

#endif // STLMesh_h
