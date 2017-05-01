// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/InitialAMR.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Initial (before t=0) adaptive mesh refinement options
  \details   Initial (before t=0) adaptive mesh refinement options
*/
// *****************************************************************************
#ifndef InciterInitialAMROptions_h
#define InciterInitialAMROptions_h

#include "NoWarning/vector.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! Mesh partitioning algorithm types
//! \author J. Bakosi
enum class InitialAMRType : uint8_t { NONE
                                    , UNIFORM };

//! \brief Pack/Unpack InitialAMRType: forward overload to generic
//!   enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, InitialAMRType& e )
{ PUP::pup( p, e ); }

//! InitialAMR options: outsource searches to base templated on enum type
//! \author J. Bakosi
class InitialAMR : public tk::Toggle< InitialAMRType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::amr_uniform >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit InitialAMR() :
      tk::Toggle< InitialAMRType >(
        //! Group, i.e., options, name
        "Initial AMR type",
        //! Enums -> names
        { { InitialAMRType::NONE, "none" },
          { InitialAMRType::UNIFORM, kw::amr_uniform::name() } },
        //! keywords -> Enums
        { { "none", InitialAMRType::NONE },
          { kw::amr_uniform::string(), InitialAMRType::UNIFORM } } ) {}
};

} // ctr::
} // tk::

#endif // InciterInitialAMROptions_h
