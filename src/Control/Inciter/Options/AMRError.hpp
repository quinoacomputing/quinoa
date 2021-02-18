// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/AMRError.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Options for computing error estimates for adaptive mesh refinement
  \details   Options for computing error estimates for adaptive mesh refinement.
*/
// *****************************************************************************
#ifndef InciterAMRErrorOptions_h
#define InciterAMRErrorOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

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
    using keywords = brigand::list< kw::amr_jump
                                  , kw::amr_hessian >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit AMRError() :
      tk::Toggle< AMRErrorType >(
        //! Group, i.e., options, name
        kw::amr_error::name(),
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
