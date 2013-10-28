//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.h
  \author    J. Bakosi
  \date      Mon Oct 28 07:21:56 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************
#ifndef Mix_h
#define Mix_h

#include <Model.h>

namespace quinoa {

//! Mix model base
class Mix : public Model {

  public:
    //! Constructor
    explicit Mix() {}
//     explicit Mix(tk::real* const particles) :
//       Model(particles,
//             base.control.get<ctr::component, ctr::npar>(),
//             base.control.nprop()),
//       m_offset(base.control.scalarOffset()),
//       m_nscalar(base.control.get<ctr::component, ctr::nscalar>()) {
//       ErrChk(m_nscalar > 0, tk::ExceptType::FATAL, "Wrong number of scalars");
//     }

    //! Destructor
    ~Mix() noexcept override = default;

//     //! Initialize particles
//     virtual void init(int p, int tid) = 0;
// 
//     //! Advance particles in mix model
//     virtual void advance(int p, int tid, tk::real dt) = 0;

  protected:
//     const int m_offset;             //!< Scalar-offset relative to base
//     const int m_nscalar;            //!< Number of mixing scalars

//     //! Initialize scalars with zero
//     void initZero(int p) {
//       memset(m_particles + p*m_nprop + m_offset, 0, m_nscalar*sizeof(tk::real));
//     }

  private:
    //! Don't permit copy constructor
    Mix(const Mix&) = delete;
    //! Don't permit copy assigment
    Mix& operator=(const Mix&) = delete;
    //! Don't permit move constructor
    Mix(Mix&&) = delete;
    //! Don't permit move assigment
    Mix& operator=(Mix&&) = delete;
};

} // quinoa::

#endif // Mix_h
