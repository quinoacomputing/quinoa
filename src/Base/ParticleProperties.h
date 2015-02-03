//******************************************************************************
/*!
  \file      src/Base/ParticleProperties.h
  \author    J. Bakosi
  \date      Sun 01 Feb 2015 07:29:59 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ParticleProperties for storing particle data.
  \details   ParticleProperties for storing particle data with various accessors
    for particle-, and property-major data layout policies. See also data layout
    [design](layout.html).
*/
//******************************************************************************
#ifndef ParticleProperties_h
#define ParticleProperties_h

#include <make_unique.h>

#include <Types.h>
#include <Config.h>

namespace tk {

//! Tags for selecting data layout policies
const uint8_t ParEqComp = 0;
const uint8_t EqCompPar = 1;

//! Zero-runtime-cost data-layout wrappers for particle properties with
//! type-based compile-time dispatch
template< uint8_t Layout >
class ParticleProperties {

  public:
    //! Constructor
    //! \param[in] npar Number of particles to allocate memory for
    //! \param[in] nprop Number properties, i.e., scalar variables, per particle
    //! \author J. Bakosi
    explicit ParticleProperties( uint64_t npar = 0, uint32_t nprop = 0) :
      m_ptr( tk::make_unique< tk::real[] >( npar*nprop ) ),
      m_npar( npar ),
      m_nprop( nprop ) {}

    //! Data access dispatch. Public interface to particle data access. Use it
    //! as ParProps( p, c, o ), where p is the particle index, c is component
    //! index specifying the scalar equation within a system of equations, and o
    //! is the offset specifying the position at which the system resides among
    //! other systems.
    //! \param[in] particle Particle index
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] offset System offset specifying the position of the system of
    //!   equations among other systems
    //! \return Reference to particle data of type tk::real
    //! \author J. Bakosi
    inline tk::real&
    operator()( int particle, int component, int offset ) const noexcept
    { return access( particle, component, offset, int2type< Layout >() ); }

    //! Const ptr to physical variable access dispatch. Public interface to
    //! physical variable access. cptr() and cvar() are intended to be used
    //! together in case component and offset would be expensive to compute for
    //! data access via the function call operator. In essence, cptr() returns
    //! part of the address known based on component and offset and intended to
    //! be used in a setup phase. Then cvar() takes this partial address and
    //! finishes the address calculation given the particle id. The following
    //! two data accesses are equivalent (modulo constness):
    //!  * real& value = operator()( par, comp, offs );
    //!  * const real* p = cptr( comp, offs );
    //!    const real& value = cvar( p, par );
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] offset System offset specifying the position of the system of
    //!   equations among other systems
    //! \return Pointer to particle data of type tk::real for use with cvar()
    //! \see Check out Statistics::setupOrdinary() and
    //!   Statistics::accumulateOrd() in Statistics/Statistics.C to see cptr()
    //!   and cvar() in action.
    //! \author J. Bakosi
    inline const tk::real*
    cptr( int component, int offset ) const noexcept
    { return cptr( component, offset, int2type< Layout >() ); }

    //! Const physical variable access dispatch. Public interface to physical
    //! variable access. cptr() and cvar() are intended to be used together in
    //! case component and offset would be expensive to compute for data access
    //! via the function call operator. In essence, cptr() returns part of the
    //! address known based on component and offset and intended to be used in a
    //! setup phase. Then cvar() takes this partial address and finishes the
    //! address calculation given the particle id. The following two data
    //! accesses are equivalent (modulo constness):
    //!  * real& value = operator()( par, comp, offs );
    //!  * const real* p = cptr( comp, offs );
    //!    const real& value = cvar( p, par );
    //! \param[in] pt Pointer to particle data of type tk::real as returned from
    //!   cptr()
    //! \param[in] particle Particle index
    //! \return Reference to particle data of type tk::real
    //! \see Check out Statistics::setupOrdinary() and
    //!   Statistics::accumulateOrd() in Statistics/Statistics.C to see cptr()
    //!   and cvar() in action.
    //! \author J. Bakosi
    inline const tk::real&
    cvar( const tk::real* const pt, int particle ) const noexcept
    { return cvar( pt, particle, int2type< Layout >() ); }

    //! Raw pointer access to particle data.
    //! \return Raw pointer to array of type tk::real
    //! \author J. Bakosi
    inline tk::real* ptr() const noexcept { return m_ptr.get(); }

    //! Total Size access.
    //! \return Total number of real numbers stored in the entire array: number
    //! of particles * number of properties/particle
    //! \author J. Bakosi
    inline uint64_t size() const noexcept { return m_npar * m_nprop; }

    //! Number of particles access.
    //! \return Number of particles
    //! \author J. Bakosi
    inline uint64_t npar() const noexcept { return m_npar; }

    //! Layout name dispatch.
    //! \return The name of data layout used
    //! \author J. Bakosi
    inline const char* major() const noexcept
    { return major( int2type< Layout >() ); }

  private:
   //! Transform a compile-time uint8_t into a type, used for dispatch
   //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
   //!   Patterns Applied, Addison-Wesley Professional, 2001.
   //! \author J. Bakosi
   template< uint8_t m > struct int2type { enum { value = m }; };

   //! Overloads for the various data accesses
   //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
   //!   Patterns Applied, Addison-Wesley Professional, 2001.
   //! \author J. Bakosi
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
   //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
   //!   Patterns Applied, Addison-Wesley Professional, 2001.
   //! \author J. Bakosi
   inline const tk::real*
   cptr( int component, int offset, int2type< ParEqComp > ) const noexcept {
     return m_ptr.get() + component + offset;
   }
   inline const tk::real*
   cptr( int component, int offset, int2type< EqCompPar > ) const noexcept {
     return m_ptr.get() + (offset+component)*m_npar;
   }

   // Overloads for the various const physical variable accesses
   //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
   //!   Patterns Applied, Addison-Wesley Professional, 2001.
   //! \author J. Bakosi
   inline const tk::real&
   cvar( const tk::real* const pt, int particle, int2type< ParEqComp > ) const
   noexcept {
     return *(pt + particle*m_nprop);
   }
   inline const tk::real&
   cvar( const tk::real* const pt, int particle, int2type< EqCompPar > ) const
   noexcept {
     return *(pt + particle);
   }

   // Overloads for the name-queries of data lauouts
   //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
   //!   Patterns Applied, Addison-Wesley Professional, 2001.
   //! \author J. Bakosi
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

//! Select data layout policy for particle properties at compile-time
#if   defined LAYOUT_PARTICLE_MAJOR
using ParProps = ParticleProperties< ParEqComp >;
#elif defined LAYOUT_EQUATION_MAJOR
using ParProps = ParticleProperties< EqCompPar >;
#endif

} // tk::

#endif // ParticleProperties_h
