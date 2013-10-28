//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Mon Oct 28 07:22:29 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <Types.h>
#include <Exception.h>

namespace quinoa {

//! Model base
class Model {

  protected:
    //! Constructor: protected, designed to be base-only
    explicit Model() {}
//     explicit Model(tk::real* const particles,
//                    const uint64_t npar,
//                    int nprop) : m_particles(particles),
//                                 m_npar(npar),
//                                 m_nprop(nprop),
//                                 m_str(nullptr) {
//       Assert(m_particles != nullptr, tk::ExceptType::FATAL,
//              "Particles pointer null?");
//       ErrChk(m_npar > 0, tk::ExceptType::FATAL,
//              "Wrong number of particles");
//       Assert(m_nprop != 0, tk::ExceptType::FATAL,
//              "Number of particle properties zero?");
//     }

    //! Destructor: protected, designed to be freed via children-only
    virtual ~Model() noexcept = default;

//     tk::real* const m_particles;    //!< Particles
//     const uint64_t m_npar;          //!< Number of particles
//     const int m_nprop;              //!< Number of particle properties

  private:
    //! Don't permit copy constructor
    Model(const Model&) = delete;
    //! Don't permit copy assigment
    Model& operator=(const Model&) = delete;
    //! Don't permit move constructor
    Model(Model&&) = delete;
    //! Don't permit move assigment
    Model& operator=(Model&&) = delete;
};

} // quinoa::

#endif // Model_h
