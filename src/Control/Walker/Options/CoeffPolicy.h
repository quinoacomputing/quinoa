//******************************************************************************
/*!
  \file      src/Control/Walker/Options/CoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 15 Dec 2015 09:50:48 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Differential equation coefficients policy options
  \details   Differential equation coefficients policy options
*/
//******************************************************************************
#ifndef CoeffPolicyOptions_h
#define CoeffPolicyOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace walker {
namespace ctr {

//! Differential equation coefficients policy types
//! \author J. Bakosi
enum class CoeffPolicyType : uint8_t { CONSTANT=0
                                     , DECAY
                                     , HOMOGENEOUS_DECAY
                                     , MONTE_CARLO_HOMOGENEOUS_DECAY
                                     , HYDROTIMESCALE_HOMOGENEOUS_DECAY
                                     };

//! Pack/Unpack CoeffPolicyType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, CoeffPolicyType& e ) { PUP::pup( p, e ); }

//! \brief CoeffPolicy options: outsource searches to base templated on enum
//!   type
//! \author J. Bakosi
class CoeffPolicy : public tk::Toggle< CoeffPolicyType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::constant
                                       , kw::decay
                                       , kw::homdecay
                                       , kw::montecarlo_homdecay
                                       , kw::hydrotimescale
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit CoeffPolicy() :
      Toggle< CoeffPolicyType >(
        //! Group, i.e., options, name
        "Coefficients Policy",
        //! Enums -> names
        { { CoeffPolicyType::CONSTANT, kw::constant::name() },
          { CoeffPolicyType::DECAY, kw::decay::name() },
          { CoeffPolicyType::HOMOGENEOUS_DECAY, kw::homdecay::name() },
          { CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY,
            kw::montecarlo_homdecay::name() },
          { CoeffPolicyType::HYDROTIMESCALE_HOMOGENEOUS_DECAY,
            kw::hydrotimescale::name() } },
        //! keywords -> Enums
        {  { kw::constant::string(), CoeffPolicyType::CONSTANT },
           { kw::decay::string(), CoeffPolicyType::DECAY },
           { kw::homdecay::string(), CoeffPolicyType::HOMOGENEOUS_DECAY },
           { kw::montecarlo_homdecay::string(),
             CoeffPolicyType::MONTE_CARLO_HOMOGENEOUS_DECAY },
           { kw::hydrotimescale::string(),
             CoeffPolicyType::HYDROTIMESCALE_HOMOGENEOUS_DECAY } } ) {}
};

} // ctr::
} // walker::

#endif // oeffPolicyOptions_h
