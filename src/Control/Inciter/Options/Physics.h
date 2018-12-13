// *****************************************************************************
/*!
  \file      src/Control/Inciter/Options/Physics.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics options for inciter
  \details   Physics options for inciter
*/
// *****************************************************************************
#ifndef PhysicsOptions_h
#define PhysicsOptions_h

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace inciter {
namespace ctr {

//! Physics types
enum class PhysicsType : uint8_t { ADVECTION=0,
                                   ADVDIFF,
                                   EULER,
                                   NAVIERSTOKES,
                                   MULTIMAT_VELEQ };

//! Pack/Unpack PhysicsType: forward overload to generic enum class packer
inline void operator|( PUP::er& p, PhysicsType& e ) { PUP::pup( p, e ); }

//! \brief Physics options: outsource to base templated on enum type
class Physics : public tk::Toggle< PhysicsType > {

  public:
    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::advection
                                  , kw::advdiff
                                  , kw::compflow_euler
                                  , kw::compflow_navierstokes
                                  , kw::multimat_compflow_veleq
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit Physics() :
      tk::Toggle< PhysicsType >(
        //! Group, i.e., options, name
        kw::physics::name(),
        //! Enums -> names (if defined, policy codes, if not, name)
        { { PhysicsType::ADVECTION, kw::advection::name() },
          { PhysicsType::ADVDIFF, kw::advdiff::name() },
          { PhysicsType::EULER, kw::compflow_euler::name() },
          { PhysicsType::NAVIERSTOKES, kw::compflow_navierstokes::name() },
          { PhysicsType::MULTIMAT_VELEQ, kw::multimat_compflow_veleq::name() }
        },
        //! keywords -> Enums
        { { kw::advection::string(), PhysicsType::ADVECTION },
          { kw::advdiff::string(), PhysicsType::ADVDIFF },
          { kw::compflow_euler::string(), PhysicsType::EULER },
          { kw::compflow_navierstokes::string(), PhysicsType::NAVIERSTOKES },
          { kw::multimat_compflow_veleq::string(),
            PhysicsType::MULTIMAT_VELEQ } } )
    {
       brigand::for_each< keywords >( assertPolicyCodes() );
    }

    //! \brief Return policy code based on Enum
    //! \param[in] p Enum value of the physics option requested
    //! \return Policy code of the option
    const std::string& code( PhysicsType p ) const {
      using tk::operator<<;
      auto it = policy.find( p );
      Assert( it != end(policy),
              std::string("Cannot find policy code for physics \"") << p <<
                "\"" );
      return it->second;
    }

  private:
    //! Function object for ensuring the existence of policy codes
    struct assertPolicyCodes {
      //! \brief Function call operator templated on the type to assert the
      //!   existence of a policy code
      template< typename U > void operator()( brigand::type_<U> ) {
        static_assert( tk::HasTypedefCode< typename U::info >::value,
                       "Policy code undefined for keyword" );
      }
    };

    //! Enums -> policy code
    std::map< PhysicsType, std::string > policy {
        { PhysicsType::ADVECTION, *kw::advection::code() }
      , { PhysicsType::ADVDIFF, *kw::advdiff::code() }
      , { PhysicsType::EULER, *kw::compflow_euler::code() }
      , { PhysicsType::NAVIERSTOKES, *kw::compflow_navierstokes::code() }
      , { PhysicsType::MULTIMAT_VELEQ, *kw::multimat_compflow_veleq::code() }
    };
};

} // ctr::
} // inciter::

#endif // PhysicsOptions_h
