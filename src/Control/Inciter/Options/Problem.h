// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Problem.h
  \author    J. Bakosi
  \date      Wed 06 Jul 2016 12:25:08 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem options for inciter
  \details   Problem options for inciter
*/
// *****************************************************************************
#ifndef ProblemOptions_h
#define ProblemOptions_h

#include <boost/mpl/vector.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Problem types
//! \author J. Bakosi
enum class ProblemType : uint8_t { USER_DEFINED=0,
                                   SHEAR_DIFF,
                                   DIR_NEU,
                                   SLOT_CYL };

//! Pack/Unpack ProblemType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, ProblemType& e ) { PUP::pup( p, e ); }

//! \brief Problem options: outsource to base templated on enum type
//! \author J. Bakosi
class Problem : public tk::Toggle< ProblemType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::shear_diff
                                       , kw::dir_neu
                                       , kw::slot_cyl
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit Problem() :
      Toggle< ProblemType >(
        //! Group, i.e., options, name
        "Test problem",
        //! Enums -> names
        { { ProblemType::USER_DEFINED, kw::user_defined::name() },
          { ProblemType::SHEAR_DIFF, kw::shear_diff::name() },
          { ProblemType::DIR_NEU, kw::dir_neu::name() },
          { ProblemType::SLOT_CYL, kw::slot_cyl::name() } },
        //! keywords -> Enums
        { { kw::user_defined::string(), ProblemType::USER_DEFINED },
          { kw::shear_diff::string(), ProblemType::SHEAR_DIFF },
          { kw::dir_neu::string(), ProblemType::DIR_NEU },
          { kw::slot_cyl::string(), ProblemType::SLOT_CYL } } ) {}
};

} // ctr::
} // inciter::

#endif // ProblemOptions_h
