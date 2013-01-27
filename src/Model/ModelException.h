//******************************************************************************
/*!
  \file      src/Model/ModelException.h
  \author    J. Bakosi
  \date      Sun 27 Jan 2013 12:18:19 PM MST
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
enum ModelExceptType { MIX_EXCEPT=0,               //!< Mix model exception
                       HYDROMODEL_EXCEPT,          //!< HydroModel exception
                       NO_SUCH_MODEL,              //!< No such model
                       ALREADY_ALLOCATED,          //!< Entry alread allocated
                       NUM_MODEL_EXCEPT
};

//! Model exception error messages
const string ModelMsg[NUM_MODEL_EXCEPT] = {
  "MixModel: ",
  "HydroModel: ",
  "No such model",
  "Memory entry already allocated"
};

//! ModelException : Exception
class ModelException : public Exception {

  public:
    //! Constructor
    ModelException(ExceptType except,
                   ModelExceptType modelExcept,
                   const string& file,
                   const string& func,
                   const unsigned int& line) :
      Exception(except, file, func, line), m_except(modelExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    ModelException(ModelException&&) = default;

    //! Don't permit copy constructor
    // ICC: should be deleted and private
    ModelException(const ModelException&);

    //! Destructor
    //virtual ~ModelException() {}

    //! Handle ModelException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy assignment
    ModelException& operator=(const ModelException&) = delete;
    //! Don't permit move assignment
    ModelException& operator=(ModelException&&) = delete;

    //! Model exception type (MIXMODEL, etc.)
    ModelExceptType m_except;
};

} // namespace Quinoa

#endif // ModelException_h
