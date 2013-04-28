//******************************************************************************
/*!
  \file      src/Model/ModelException.h
  \author    J. Bakosi
  \date      Sat 27 Apr 2013 08:28:02 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ModelException
  \details   ModelException
*/
//******************************************************************************
#ifndef ModelException_h
#define ModelException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! Model exception types
enum ModelExceptType { MIX_EXCEPT=0,           //!< Mix model exception
                       HYDROMODEL_EXCEPT,      //!< HydroModel exception
                       NO_SUCH_MODEL,          //!< No such model
                       BAD_NPAR,               //!< Wrong number of particles
                       ALREADY_ALLOCATED,      //!< Entry alread allocated
                       NUM_MODEL_EXCEPT
};

//! Model exception error messages
const string ModelMsg[NUM_MODEL_EXCEPT] = {
  "MixModel: ",
  "HydroModel: ",
  "No such model",
  "Wrong number of particles",
  "Memory entry already allocated"
};

//! ModelException : Exception
class ModelException : public Exception {

  public:
    //! Constructor
    explicit ModelException(const ExceptType except,
                            const ModelExceptType modelExcept,
                            const string& file,
                            const string& func,
                            const unsigned int& line,
                            const string& message = "") noexcept :
      Exception(except,
                file,
                func,
                line,
                ModelMsg[static_cast<int>(modelExcept)] + message) {}

    //! Destructor
    virtual ~ModelException() noexcept = default;

    //! Move constructor for throws, default compiler generated
    ModelException(ModelException&&) = default;

  protected:
    //! Permit copy constructor only for children
    ModelException(const ModelException&) = default;

  private:
    //! Don't permit copy assignment
    ModelException& operator=(const ModelException&) = delete;
    //! Don't permit move assignment
    ModelException& operator=(ModelException&&) = delete;
};

} // namespace Quinoa

#endif // ModelException_h
