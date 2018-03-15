// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/AMRError.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Options for computing error estimates for adaptive mesh refinement
  \details   Options for computing error estimates for adaptive mesh refinement.
*/
// *****************************************************************************
#ifndef InciterAMRErrorOptions_h
#define InciterAMRErrorOptions_h

#include "NoWarning/vector.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Mesh partitioning algorithm types
enum class AMRErrorType : uint8_t { JUMP
                                  , HESSIAN };

//! Pack/Unpack AMRErrorType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, AMRErrorType& e ) { PUP::pup( p, e ); }

//! AMRError options: outsource searches to base templated on enum type
class AMRError : public tk::Toggle< AMRErrorType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::amr_jump
                                       , kw::amr_hessian >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit AMRError() :
      tk::Toggle< AMRErrorType >(
        //! Group, i.e., options, name
        "Error estimator",
        //! Enums -> names
        { { AMRErrorType::JUMP, kw::amr_jump::name() },
          { AMRErrorType::HESSIAN, kw::amr_hessian::name() } },
        //! keywords -> Enums
        { { kw::amr_jump::string(), AMRErrorType::JUMP },
          { kw::amr_hessian::string(), AMRErrorType::HESSIAN } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterAMRErrorOptions_h
