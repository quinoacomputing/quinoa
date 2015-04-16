//******************************************************************************
/*!
  \file      src/Control/Inciter/Options/Problem.h
  \author    J. Bakosi
  \date      Wed 15 Apr 2015 10:45:05 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Test problem options for inciter
  \details   Test problem options for inciter
*/
//******************************************************************************
#ifndef ProblemOptions_h
#define ProblemOptions_h

#include <boost/mpl/vector.hpp>

#include <Toggle.h>
#include <Keywords.h>
#include <PUPUtil.h>

namespace inciter {
namespace ctr {

//! Test problem types
//! \author J. Bakosi
enum class ProblemType : uint8_t { NO_PROB=0,
                                   SHEAR_DIFF,
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
        { { ProblemType::SHEAR_DIFF, kw::shear_diff::name() },
          { ProblemType::SLOT_CYL, kw::slot_cyl::name() } },
        //! keywords -> Enums
        { { kw::shear_diff::string(), ProblemType::SHEAR_DIFF },
          { kw::slot_cyl::string(), ProblemType::SLOT_CYL } } ) {}
};

} // ctr::
} // inciter::

#endif // ProblemOptions_h
