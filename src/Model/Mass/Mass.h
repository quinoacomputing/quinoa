//******************************************************************************
/*!
  \file      src/Model/Mass/Mass.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:15:11 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mass model base
  \details   Mass mode lbase
*/
//******************************************************************************
#ifndef Mass_h
#define Mass_h

#include <cstring>

#include <QuinoaTypes.h>
#include <Model.h>
#include <QuinoaControl.h>

namespace quinoa {

//! Mass model base
class Mass : public Model {

  public:
    //! Constructor
    explicit Mass(const Base& base, real* const particles) :
      Model(base,
            particles,
            base.control.get<ctr::component, ctr::npar>(),
            base.control.nprop()),
      m_offset(0),
      m_ndensity(base.control.get<ctr::component, ctr::ndensity>()) {
      ErrChk(m_ndensity > 0, ExceptType::FATAL,
             "Wrong number of particle densities");
    }

    //! Destructor
    ~Mass() noexcept override = default;

    //! Initialize particles
    virtual void init() = 0;

    //! Advance particles in mass model
    virtual void advance(const real& dt) = 0;

  protected:
    const int m_offset;             //!< Mass-offset relative to base
    const int m_ndensity;           //!< Number of density components

//     //! Initialize densities with beta symmetric PDF
//     //! \param[in] alpha  First shape parameter
//     //! \param[in] beta   Second shape parameter
//     //! \param[in] disp   Displacement (i.e., shift) parameter
//     //! \param[in] scale  Scale parameter
//     void initBeta(const real alpha,
//                   const real beta,
//                   const real disp,
//                   const real scale) {
//       for (uint64_t p=0; p<m_npar; ++p) {
//         m_rndStr->beta(VSL_RNG_METHOD_BETA_CJA,
//                        m_str[0], 1, m_densities+p, alpha, beta, disp, scale);
//       }
//     }

  private:
    //! Don't permit copy constructor
    Mass(const Mass&) = delete;
    //! Don't permit copy assigment
    Mass& operator=(const Mass&) = delete;
    //! Don't permit move constructor
    Mass(Mass&&) = delete;
    //! Don't permit move assigment
    Mass& operator=(Mass&&) = delete;
};

} // namespace quinoa

#endif // Mass_h
