//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 12:37:13 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <Control.h>

namespace Quinoa {

class MixModel;

//! Model base
class Model {

  public:
    //! Constructor
    Model(const ModelType model, const int npel);

    //! Destructor
    ~Model();

    //! Echo informaion on model
    void echo();

    //! Initialize model
    void init();

    //! Constant accessor for number of particles/element
    const int& npel() { return m_npel; }

  private:
    //! Don't permit copy constructor
    Model(const Model&) = delete;
    //! Don't permit copy assigment
    Model& operator=(const Model&) = delete;
    //! Don't permit move constructor
    Model(Model&&) = delete;
    //! Don't permit move assigment
    Model& operator=(Model&&) = delete;

    const ModelType m_model;          //!< Model type
    const int m_npel;                 //!< Number of particles/element

    MixModel* m_mixModel;             //!< Pointer to MixModel object
};

} // namespace Quinoa

#endif // Model_h
