// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/AMRInitial.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Initial (before t=0) adaptive mesh refinement (AMR) options
  \details   Initial (before t=0) adaptive mesh refinement (AMR) options
*/
// *****************************************************************************
#ifndef InciterAMRInitialOptions_h
#define InciterAMRInitialOptions_h

#include "NoWarning/vector.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Initial AMR types
enum class AMRInitialType : uint8_t { NONE
                                    , UNIFORM };

//! Pack/Unpack AMRInitialType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, AMRInitialType& e )
{ PUP::pup( p, e ); }

//! AMRInitial options: outsource searches to base templated on enum type
class AMRInitial : public tk::Toggle< AMRInitialType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = boost::mpl::vector< kw::amr_uniform >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit AMRInitial() :
      tk::Toggle< AMRInitialType >(
        //! Group, i.e., options, name
        "Initial AMR type",
        //! Enums -> names
        { { AMRInitialType::NONE, "none" },
          { AMRInitialType::UNIFORM, kw::amr_uniform::name() } },
        //! keywords -> Enums
        { { "none", AMRInitialType::NONE },
          { kw::amr_uniform::string(), AMRInitialType::UNIFORM } } ) {}
};

} // ctr::
} // inciter::

#endif // InciterAMRInitialOptions_h
