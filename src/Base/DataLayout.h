//******************************************************************************
/*!
  \file      src/Base/DataLayout.h
  \author    J. Bakosi
  \date      Sun 31 Jan 2016 07:16:42 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Generic data access abstraction for different data layouts
  \details   Generic data access abstraction for different data layouts. See
    also Base/DataLayout.h and the rationale discussed in the
    [design](layout.html) document.
*/
//******************************************************************************
#ifndef DataLayout_h
#define DataLayout_h

#include "Types.h"
#include "Make_unique.h"
#include "Keywords.h"

namespace tk {

//! Tags for selecting data layout policies
const uint8_t UnkEqComp = 0;
const uint8_t EqCompUnk = 1;

//! Zero-runtime-cost data-layout wrappers with type-based compile-time dispatch
template< uint8_t Layout >
class DataLayout {

  private:
    //! \brief Inherit type of number of components from keyword 'ncomp', used
    //!    also for type of offset
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor
    //! \param[in] nunk Number of unknowns to allocate memory for
    //! \param[in] nprop Number properties, i.e., scalar variables, per unknown
    //! \author J. Bakosi
    explicit DataLayout( ncomp_t nunk = 0, ncomp_t nprop = 0) :
      m_ptr( tk::make_unique< tk::real[] >( nunk*nprop ) ),
      m_nunk( nunk ),
      m_nprop( nprop ) {}

    //! \brief Data access dispatch.
    //! \details Public interface to data access. Use it as DataLayout(p,c,o),
    //!   where p is the unknown index, c is component index specifying the
    //!   scalar equation within a system of equations, and o is the offset
    //!   specifying the position at which the system resides among other
    //!   systems.
    //! \param[in] unknown Unknown index
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] offset System offset specifying the position of the system of
    //!   equations among other systems
    //! \return Reference to data of type tk::real
    //! \author J. Bakosi
    inline tk::real&
    operator()( ncomp_t unknown, ncomp_t component, ncomp_t offset )
    const noexcept {
      return access( unknown, component, offset, int2type< Layout >() );
    }

    //! \brief Const ptr to physical variable access dispatch.
    //! \details Public interface to physical variable access. cptr() and
    //!   cvar() are intended to be used together in case component and offset
    //!   would be expensive to compute for data access via the function call
    //!   operator. In essence, cptr() returns part of the address known based
    //!   on component and offset and intended to be used in a setup phase. Then
    //!   cvar() takes this partial address and finishes the address calculation
    //!   given the unknown id. The following two data accesses are equivalent
    //!   (modulo constness):
    //!   * real& value = operator()( unk, comp, offs );
    //!   * const real* p = cptr( comp, offs );
    //!     const real& value = cvar( p, unk );
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] offset System offset specifying the position of the system of
    //!   equations among other systems
    //! \return Pointer to data of type tk::real for use with cvar()
    //! \see Client code for cptr() and cvar() in Statistics::setupOrdinary()
    //!   and Statistics::accumulateOrd() in Statistics/Statistics.C.
    //! \author J. Bakosi
    inline const tk::real*
    cptr( ncomp_t component, ncomp_t offset ) const noexcept
    { return cptr( component, offset, int2type< Layout >() ); }

    //! \brief Const physical variable access dispatch.
    //! \details Public interface to physical variable access. cptr() and cvar()
    //!   are intended to be used together in case component and offset would be
    //!   expensive to compute for data access via the function call operator.
    //!   In essence, cptr() returns part of the address known based on
    //!   component and offset and intended to be used in a setup phase. Then
    //!   cvar() takes this partial address and finishes the address calculation
    //!   given the unknown id. The following two data accesses are equivalent
    //!   (modulo constness):
    //!   * real& value = operator()( unk, comp, offs );
    //!   * const real* p = cptr( comp, offs );
    //!    const real& value = cvar( p, unk );
    //! \param[in] pt Pointer to data of type tk::real as returned from cptr()
    //! \param[in] unknown Unknown index
    //! \return Reference to data of type tk::real
    //! \see Client code for cptr() and cvar() in Statistics::setupOrdinary()
    //!   and Statistics::accumulateOrd() in Statistics/Statistics.C.
    //! \author J. Bakosi
    inline const tk::real&
    cvar( const tk::real* const pt, ncomp_t unknown ) const noexcept
    { return cvar( pt, unknown, int2type< Layout >() ); }

    //! Raw pointer access to data.
    //! \return Raw pointer to array of type tk::real
    //! \author J. Bakosi
    inline tk::real* ptr() const noexcept { return m_ptr.get(); }

    //! Total Size access.
    //! \return Total number of real numbers stored in the entire array: number
    //! of unknowns * number of properties/unknown.
    //! \author J. Bakosi
    inline ncomp_t size() const noexcept { return m_nunk * m_nprop; }

    //! Number of unknown access.
    //! \return Number of unknowns
    //! \author J. Bakosi
    inline ncomp_t nunk() const noexcept { return m_nunk; }

    //! Layout name dispatch.
    //! \return The name of the data layout used
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
    //! \param[in] unknown Unknown index
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] offset System offset specifying the position of the system of
    //!   equations among other systems
    //! \return Reference to data of type tk::real
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    //! \author J. Bakosi
    inline tk::real&
    access( ncomp_t unknown, ncomp_t component, ncomp_t offset,
            int2type< UnkEqComp > ) const noexcept
    {
      return *(m_ptr.get() + unknown*m_nprop + offset + component);
    }
    inline tk::real&
    access( ncomp_t unknown, ncomp_t component, ncomp_t offset,
            int2type< EqCompUnk > ) const noexcept
    {
      return *(m_ptr.get() + (offset+component)*m_nunk + unknown);
    }

    // Overloads for the various const ptr to physical variable accesses
    //! \param[in] component Component index, i.e., position of a scalar within
    //!   a system
    //! \param[in] offset System offset specifying the position of the system of
    //!   equations among other systems
    //! \return Pointer to data of type tk::real for use with cvar()
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    //! \author J. Bakosi
    inline const tk::real*
    cptr( ncomp_t component, ncomp_t offset, int2type< UnkEqComp > ) const
    noexcept {
      return m_ptr.get() + component + offset;
    }
    inline const tk::real*
    cptr( ncomp_t component, ncomp_t offset, int2type< EqCompUnk > ) const
    noexcept {
      return m_ptr.get() + (offset+component)*m_nunk;
    }

    // Overloads for the various const physical variable accesses
    //! \param[in] pt Pointer to data of type tk::real as returned from cptr()
    //! \param[in] unknown Unknown index
    //! \return Reference to data of type tk::real
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    //! \author J. Bakosi
    inline const tk::real&
    cvar( const tk::real* const pt, ncomp_t unknown, int2type< UnkEqComp > )
    const noexcept {
      return *(pt + unknown*m_nprop);
    }
    inline const tk::real&
    cvar( const tk::real* const pt, ncomp_t unknown, int2type< EqCompUnk > )
    const noexcept {
      return *(pt + unknown);
    }

    // Overloads for the name-queries of data lauouts
    //! \return The name of the data layout used
    //! \see A. Alexandrescu, Modern C++ Design: Generic Programming and Design
    //!   Patterns Applied, Addison-Wesley Professional, 2001.
    //! \author J. Bakosi
    inline const char* major( int2type< UnkEqComp > ) const noexcept {
      return "unknown-major";
    }
    inline const char* major( int2type< EqCompUnk > ) const noexcept {
      return "equation-major";
    }

    const std::unique_ptr< tk::real[] > m_ptr; //!< Particle data pointer
    const ncomp_t m_nunk;                      //!< Number of unknowns
    const ncomp_t m_nprop;                     //!< Number of properties/unknown
};

} // tk::

#endif // DataLayout_h
