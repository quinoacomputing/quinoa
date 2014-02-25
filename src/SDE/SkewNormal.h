//******************************************************************************
/*!
  \file      src/SDE/SkewNormal.h
  \author    J. Bakosi
  \date      Mon 24 Feb 2014 08:12:33 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Skew-normal SDE
  \details   Skew-normal SDE
*/
//******************************************************************************
#ifndef SkewNormal_h
#define SkewNormal_h

#include <SDE.h>
#include <SkewNormalCoeffPolicy.h>

namespace quinoa {

//! SkewNormal
template< class Init, class Coefficients >
class SkewNormal : public SDE< Init, Coefficients > {

  public:
    //! SDE base shorthand
    using sde = SDE< Init >;

    //! Constructor
    explicit SkewNormal( const Base& base,
                         const ParProps& particles,
                         int offset,
                         int ncomp ) :
      sde( base,
           base.control.get< tag::param, tag::dirichlet, tk::tag::rng >(),
           particles,
           offset,
           ncomp ),
      m_coeff( m_coeffPolicy,
               base.control.get< tag::param, tag::skewnormal, tag::sigma >(),
               base.control.get< tag::param, tag::skewnormal, tag::timescale >(),
               base.control.get< tag::param, tag::lambda, tag::lambda >(),
               m_sigma, m_timescale, m_lambda ) {}

    //! Return coefficients policy
    const std::string& coeffPolicy() const noexcept override {
      return m_coeffPolicy;
    }

    //! Pull base class data to scope
    using sde::m_particles;
    using sde::m_offset;
    using sde::m_ncomp;
    using sde::m_rng;

    //! Advance particles
    void advance(int p, int tid, tk::real dt) override {
    }

  private:
    //! Don't permit copy constructor
    SkewNormal(const SkewNormal&) = delete;
    //! Don't permit copy assigment
    SkewNormal& operator=(const SkewNormal&) = delete;
    //! Don't permit move constructor
    SkewNormal(SkewNormal&&) = delete;
    //! Don't permit move assigment
    SkewNormal& operator=(SkewNormal&&) = delete;

    tk::real m_sigma;                   //!< SDE coefficients
    tk::real m_timescale;
    tk::real m_lambda;

    std::string m_coeffPolicy;          //!< Coefficients policy name
    Coefficients m_coeff;               //!< Coefficients policy
};

} // quinoa::

#endif // SkewNormal_h
