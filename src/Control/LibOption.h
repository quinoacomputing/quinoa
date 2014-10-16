//******************************************************************************
/*!
  \file      src/Control/LibOption.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 10:04:45 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Library toggle
  \details   Library toggle
*/
//******************************************************************************
#ifndef LibOption_h
#define LibOption_h

namespace tk {

//! Library toggle
template< class Option >
class LibOption {

  public:
    //! Lookup library option parameter
    const typename Option::ParamType& param( typename Option::EnumType value )
    const { return Option().param( value ); }

    //! Lookup library option library
    typename Option::LibType lib( typename Option::EnumType value ) const
    { return Option().lib( value ); }
};

} // tk::

#endif // LibOption_h
