// *****************************************************************************
/*!
  \file      src/Control/Walker/Options/VelocityVariant.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Velocity model variants
  \details   Velocity model variants for walker.
*/
// *****************************************************************************
#ifndef VelocityVariantOptions_h
#define VelocityVariantOptions_h

#include <brigand/sequences/list.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace walker {
namespace ctr {

//! Velocity model variant types
enum class VelocityVariantType : uint8_t { SLM=0
                                         , GLM };

//! \brief Pack/Unpack VelocityVariantType: forward overload to generic enum
//!    class packer
inline void operator|( PUP::er& p, VelocityVariantType& e )
{ PUP::pup( p, e ); }

//! Velocity model variants: outsource to base templated on enum type
class VelocityVariant : public tk::Toggle< VelocityVariantType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::slm
                                  , kw::glm
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit VelocityVariant() :
      tk::Toggle< VelocityVariantType >(
        //! Group, i.e., options, name
        "Model variant",
        //! Enums -> names
        { { VelocityVariantType::SLM, kw::slm::name() },
          { VelocityVariantType::GLM, kw::glm::name() } },
        //! keywords -> Enums
        { { kw::slm::string(), VelocityVariantType::SLM },
          { kw::glm::string(), VelocityVariantType::GLM } } )
    {}
};

} // ctr::
} // walker::

#endif // VelocityVariantOptions_h
