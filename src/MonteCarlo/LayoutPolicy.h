//******************************************************************************
/*!
  \file      src/MonteCarlo/LayoutPolicy.h
  \author    J. Bakosi
  \date      Mon 27 Jan 2014 02:14:26 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Particle-, and property-major data layout policies
  \details   Particle-, and property-major data layout policies
*/
//******************************************************************************
#ifndef LayoutPolicy_h
#define LayoutPolicy_h

#include<make_unique.h>

namespace quinoa {

//! Tags for selecting particle-, or property-major data layout policies
const uint8_t ParEqComp = 0;
const uint8_t ParCompEq = 1;
const uint8_t EqCompPar = 2;
const uint8_t EqParComp = 3;
const uint8_t CompEqPar = 4;
const uint8_t CompParEq = 5;

//! Zero-runtime-cost data-layout wrappers for particle properties with
//! type-based compile-time dispatch
template< uint8_t Layout >
class ParticleProperties {

  private:
    //! Don't permit copy constructor
    ParticleProperties(const ParticleProperties&) = delete;
    //! Don't permit copy assigment
    ParticleProperties& operator=(const ParticleProperties&) = delete;
    //! Don't permit move constructor
    ParticleProperties(ParticleProperties&&) = delete;
    //! Don't permit move assigment
    ParticleProperties& operator=(ParticleProperties&&) = delete;

   //! Transform a compile-time uint8_t into a type
   template< uint8_t m >
   struct int2type {
     enum { value = m };
   };

   // Overloads for the various data accesses
   inline tk::real& access( int particle, int property, int offset,
                            int2type< ParEqComp > ) const noexcept {
     return *(m_ptr.get() + particle*m_nprop + offset + property);
   }
   inline tk::real& access( int particle, int property, int offset,
                            int2type< ParCompEq > ) const noexcept {
     return *(m_ptr.get() + particle*m_nprop + offset + property);
   }

   // Overloads for the various const ptr accesses
   inline const tk::real*
   cptr_access( int property, int offset, int2type< ParEqComp > ) const
   noexcept {
     return m_ptr.get() + offset + property;
   }
   inline const tk::real*
   cptr_access( int property, int offset, int2type< ParCompEq > ) const
   noexcept {
     return m_ptr.get() + offset + property;
   }

   // Overloads for the name-queries of data lauouts
   inline const char* major( int2type< ParEqComp > ) const noexcept {
     return "[ particle ] [ equation ] [ component ]";
   }
   inline const char* major( int2type< ParCompEq > ) const noexcept {
     return "[ particle ] [ component ] [ equation ]";
   }
   inline const char* major( int2type< EqCompPar > ) const noexcept {
     return "[ equation ] [ component ] [ particle ]";
   }
   inline const char* major( int2type< EqParComp > ) const noexcept {
     return "[ equation ] [ particle ] [ component ]";
   }
   inline const char* major( int2type< CompEqPar > ) const noexcept {
     return "[ component ] [ equation ] [ particle ]";
   }
   inline const char* major( int2type< CompParEq > ) const noexcept {
     return "[ component ] [ particle ] [ equation ]";
   }

   const std::unique_ptr< tk::real[] > m_ptr; //!< Particle data pointer
   const int m_nprop;                         //!< Number of particle properties

  public:
    //! Constructor
    ParticleProperties( uint64_t npar, int nprop ) :
      m_ptr( tk::make_unique< tk::real[] >( npar*nprop ) ),
      m_nprop( nprop ) {}

    //! Data access dispatch
    inline tk::real&
    operator()( int particle, int property, int offset ) const noexcept {
      return access( particle, property, offset, int2type< Layout >() );
    }

    //! Const ptr access dispatch
    inline const tk::real*
    cptr( int property, int offset ) const noexcept {
      return cptr_access( property, offset, int2type< Layout >() );
    }

    //! Ptr access dispatch
    inline tk::real* ptr() const noexcept { return m_ptr.get(); }

    //! Layout name dispatch
    inline const char* major() const noexcept {
      return major( int2type< Layout >() );
    }
};

} // quinoa::

#endif // LayoutPolicy_h
