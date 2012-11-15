//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Thu Nov 15 13:27:34 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <iostream>

#include <Control.h>

using namespace std;

namespace Quinoa {

class Memory;
class MemoryEntry;
class MixModel;
class Paradigm;
class MKLRandom;

//! Model base
class Model {

  public:
    //! Constructor
    Model(const ModelType model,
          const int npel,
          Memory* memory,
          Paradigm* paradigm);

    //! Destructor
    ~Model();

    //! Echo informaion on model
    void echo();

    //! Initialize model
    void init();

    //! Allocate array for storing the element IDs of particles
    void allocNpel();

    //! Constant accessor to number of particles/element
    const int& npel() const { return m_npel; }

    //! Constant accessor to number of elements
    const int& nel() const { return m_nel; }

    //! Accessor to memory object pointer
    Memory* memory() const { return m_memory; }

    //! Accessor to parallel programming object pointer
    Paradigm* paradigm() const { return m_paradigm; }

    //! Accessor to random number generator object
    MKLRandom* getRandom() const { return m_random; }

  private:
    //! Don't permit copy constructor
    Model(const Model&) = delete;
    //! Don't permit copy assigment
    Model& operator=(const Model&) = delete;
    //! Don't permit move constructor
    Model(Model&&) = delete;
    //! Don't permit move assigment
    Model& operator=(Model&&) = delete;

    const ModelType m_model;      //!< Model type
    const int m_npel;             //!< Number of particles/element
    Memory* m_memory;             //!< Memory object pointer
    Paradigm* m_paradigm;         //!< Parallel programming object pointer
    string m_name;                //!< Name of model
    int m_nel;                    //!< Number of elements
    MKLRandom* m_random;          //!< Pointer to random number generator object
    MixModel* m_mixModel;         //!< Pointer to MixModel object

    MemoryEntry* m_elp;           //!< Array storing element ID of particle
};

} // namespace Quinoa

#endif // Model_h
