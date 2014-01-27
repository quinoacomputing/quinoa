//******************************************************************************
/*!
  \file      src/MonteCarlo/LayoutPolicy.h
  \author    J. Bakosi
  \date      Mon 27 Jan 2014 07:50:58 AM MST
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
const bool ParticleMajor = true;
const bool PropertyMajor = false;

//! Zero-runtime-cost data-layout wrappers for particle properties with
//! type-based compile-time dispatch
template< bool Major >
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

   //! Transform a compile-time bool into a type
   template< bool m >
   struct int2type {
     enum { value = m };
   };

   // Overloads for particle-, and property-major data accesses
   inline tk::real& access( int particle, int property, int offset,
                            int2type<ParticleMajor> ) const
   {
     return *(m_ptr.get() + particle*m_nprop + offset + property);
   }
   inline tk::real& access( int particle, int property, int offset,
                            int2type<PropertyMajor> ) const
   {
     return *(m_ptr.get() + particle*m_nprop + offset + property);
   }

   // Overloads for particle-, and property-major const ptr accesses
   inline const tk::real*
   cptr_access( int property, int offset, int2type<ParticleMajor> ) const {
     return m_ptr.get() + offset + property;
   }
   inline const tk::real*
   cptr_access( int property, int offset, int2type<PropertyMajor> ) const {
     return m_ptr.get() + offset + property;
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
    operator()( int particle, int property, int offset ) const {
      return access( particle, property, offset, int2type<Major>() );
    }

    //! Const ptr access dispatch
    inline const tk::real*
    cptr( int property, int offset ) const {
      return cptr_access( property, offset, int2type<Major>() );
    }

    //! Ptr access dispatch
    inline tk::real* ptr() const { return m_ptr.get(); }
};

} // quinoa::

#endif // LayoutPolicy_h
