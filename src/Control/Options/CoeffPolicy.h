//******************************************************************************
/*!
  \file      src/Control/Options/CoeffPolicy.h
  \author    J. Bakosi
  \date      Fri 23 Jan 2015 06:44:59 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Differential equation coefficients policy options
  \details   Differential equation coefficients policy options
*/
//******************************************************************************
#ifndef CoeffPolicyOptions_h
#define CoeffPolicyOptions_h

#include <boost/mpl/vector.hpp>

#include <Toggle.h>
#include <Keywords.h>
#include <PUPUtil.h>

namespace tk {
namespace ctr {

//! Differential equation coefficients policy types
//! \author J. Bakosi
enum class CoeffPolicyType : uint8_t { CONSTANT=0
                                     , JRRJ
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
                                       , kw::jrrj
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
          { CoeffPolicyType::JRRJ, kw::jrrj::name() } },
        //! keywords -> Enums
        {  { kw::constant::string(), CoeffPolicyType::CONSTANT },
           { kw::jrrj::string(), CoeffPolicyType::JRRJ } } ) {}
};

} // ctr::
} // tk::

#endif // oeffPolicyOptions_h
