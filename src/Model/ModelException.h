//******************************************************************************
/*!
  \file      src/Model/ModelException.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:00:09 2013
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
                            const unsigned int& line) :
      Exception(except, file, func, line), m_except(modelExcept) {}

    //! Destructor
    virtual ~ModelException() noexcept = default;

    //! Handle ModelException
    virtual ErrCode handleException(Driver* const driver);

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

    //! Model exception type (MIXMODEL, etc.)
    const ModelExceptType m_except;
};

} // namespace Quinoa

#endif // ModelException_h
