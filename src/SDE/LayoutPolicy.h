//******************************************************************************
/*!
  \file      src/SDE/LayoutPolicy.h
  \author    J. Bakosi
  \date      Wed Jan 15 11:29:18 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Data layout policies
  \details   Data layout policies
*/
//******************************************************************************
#ifndef LayoutPolicy_h
#define LayoutPolicy_h

namespace quinoa {

struct Data {
  Data( tk::real* const ptr, int nprop, int offset ) :
    m_ptr(ptr), m_nprop(nprop), m_offset(offset) {}

  tk::real& operator()( int particle, int property ) {
    return *(m_ptr + particle*m_nprop + m_offset);
  }

  tk::real* const m_ptr;
  const int m_nprop;
  const int m_offset;
};

//! Particle-major data layout policy
struct LayoutParticleMajor {
};

//! Property-major data layout policy
struct LayoutPropertyMajor {
};

} // quinoa::

#endif // LayoutPolicy_h
