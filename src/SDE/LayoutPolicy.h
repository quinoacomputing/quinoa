//******************************************************************************
/*!
  \file      src/SDE/LayoutPolicy.h
  \author    J. Bakosi
  \date      Tue 21 Jan 2014 10:05:40 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Data layout policies
  \details   Data layout policies
*/
//******************************************************************************
#ifndef LayoutPolicy_h
#define LayoutPolicy_h

namespace quinoa {

//! Tags to select particle-, or property-major data layout policies
const bool ParticleMajor = true;
const bool PropertyMajor = false;

template< bool Major >
class Data {

  private:
   //! Transform a compile-time bool into a type
   template< bool m >
   struct int2type {
     enum { value = m };
   };

   // Overloads for particle-, and property-major accesses
   tk::real& access( int particle, int property, int2type<ParticleMajor> ) {
     std::cout << "particle\n";
     return *(m_ptr + particle*m_nprop + m_offset + property);
   }
   tk::real& access( int particle, int property, int2type<PropertyMajor> ) {
     std::cout << "property\n";
     return *(m_ptr + particle*m_nprop + m_offset + property);
   }

   tk::real* const m_ptr;
   const int m_nprop;
   const int m_offset;

  public:
    //! Constructor
    Data( tk::real* const ptr, int nprop, int offset ) :
      m_ptr(ptr), m_nprop(nprop), m_offset(offset) {}

    //! Access dispatch
    tk::real& operator()( int particle, int property ) {
      return access( particle, property, int2type<Major>() );
    }
};

} // quinoa::

#endif // LayoutPolicy_h
