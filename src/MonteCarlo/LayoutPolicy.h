//******************************************************************************
/*!
  \file      src/MonteCarlo/LayoutPolicy.h
  \author    J. Bakosi
  \date      Tue 28 Jan 2014 02:19:06 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Particle-, and property-major data layout policies
  \details   Particle-, and property-major data layout policies
*/
//******************************************************************************
#ifndef LayoutPolicy_h
#define LayoutPolicy_h

#include<make_unique.h>

namespace quinoa {

//! Tags for selecting data layout policies
const uint8_t ParEqComp = 0;
const uint8_t EqCompPar = 1;

//! Zero-runtime-cost data-layout wrappers for particle properties with
//! type-based compile-time dispatch
template< uint8_t Layout >
class ParticleProperties {

  public:
    //! Constructor
    ParticleProperties( uint64_t npar, uint32_t nprop ) :
      m_ptr( tk::make_unique< tk::real[] >( npar*nprop ) ),
      m_npar( npar ),
      m_nprop( nprop ) {}

    //! Data access dispatch
    inline tk::real&
    operator()( int particle, int component, int offset ) const noexcept {
      return access( particle, component, offset, int2type< Layout >() );
    }

    // cptr() and cvar() are intended to be used together in case component and
    // offset would be expensive to compute for data access via the function
    // call operator. In essence, cptr() returns part of the address known
    // based on component and offset and intended to be used in a setup phase.
    // Then cvar() takes this partial address and finishes the address
    // calculation given the particle id. The following two data accesses are
    // equivalent (modulo constness):
    //  * real& value = operator()( par, comp, offs );
    //  * const real* p = cptr( comp, offs );
    //    const real& value = cvar( p, par );

    //! Const ptr to physical variable access dispatch
    inline const tk::real*
    cptr( int component, int offset ) const noexcept {
      return cptr( component, offset, int2type< Layout >() );
    }

    //! Const physical variable access dispatch
    inline const tk::real&
    cvar( const tk::real* const ptr, int particle ) const noexcept {
      return cvar( ptr, particle, int2type< Layout >() );
    }

    //! Ptr access
    inline tk::real* ptr() const noexcept { return m_ptr.get(); }

    //! Size access
    inline uint64_t size() const noexcept { return m_npar * m_nprop; }

    //! Layout name dispatch
    inline const char* major() const noexcept {
      return major( int2type< Layout >() );
    }

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
   inline tk::real&
   access( int particle, int component, int offset, int2type< ParEqComp > )
   const noexcept {
     return *(m_ptr.get() + particle*m_nprop + offset + component);
   }
   inline tk::real&
   access( int particle, int component, int offset, int2type< EqCompPar > )
   const noexcept {
     return *(m_ptr.get() + (offset+component)*m_npar + particle);
   }

   // Overloads for the various const ptr to physical variable accesses
   inline const tk::real*
   cptr( int component, int offset, int2type< ParEqComp > ) const noexcept {
     return m_ptr.get() + component + offset;
   }
   inline const tk::real*
   cptr( int component, int offset, int2type< EqCompPar > ) const noexcept {
     return m_ptr.get() + (offset+component)*m_npar;
   }

   // Overloads for the various const physical variable accesses
   inline const tk::real&
   cvar( const tk::real* const ptr, int particle, int2type< ParEqComp > ) const
   noexcept {
     return *(ptr + particle*m_nprop);
   }
   inline const tk::real&
   cvar( const tk::real* const ptr, int particle, int2type< EqCompPar > ) const
   noexcept {
     return *(ptr + particle);
   }

   // Overloads for the name-queries of data lauouts
   inline const char* major( int2type< ParEqComp > ) const noexcept {
     return "particle-major";
   }
   inline const char* major( int2type< EqCompPar > ) const noexcept {
     return "equation-major";
   }

   const std::unique_ptr< tk::real[] > m_ptr; //!< Particle data pointer
   const uint64_t m_npar;                     //!< Number of particles
   const uint32_t m_nprop;                    //!< Number of particle properties
};

} // quinoa::

#endif // LayoutPolicy_h
