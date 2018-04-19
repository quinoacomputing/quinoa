// *****************************************************************************
/*!
  \file      src/Base/Has.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     "Has-a" utilities for detecting class internals
  \details   "Has-a" utilities for detecting class internals
*/
// *****************************************************************************
#ifndef Has_h
#define Has_h

namespace tk {

//! \brief Detect if T defines type "Proxy"
//! \see http://stackoverflow.com/a/7235647
template< typename T >
struct HasTypedefProxy {
  private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;
    template<typename C> static yes test(typename C::Proxy*);
    template<typename C> static no  test(...);
  public:
    static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
};

//! \brief Detect if T defines type "alias"
//! \see http://stackoverflow.com/a/7235647
template< typename T >
struct HasTypedefAlias {
  private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;
    template<typename C> static yes test(typename C::alias*);
    template<typename C> static no  test(...);
  public:
    static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
};

//! \brief Detect if T defines type "code"
//! \see http://stackoverflow.com/a/7235647
template< typename T >
struct HasTypedefCode {
  private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;
    template<typename C> static yes test(typename C::code*);
    template<typename C> static no  test(...);
  public:
    static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
};

//! \brief Detect if T defines type "expect::type"
//! \see http://stackoverflow.com/a/7235647
template< typename T >
struct HasTypedefExpectType {
  private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;
    template<typename C> static yes test(typename C::expect::type*);
    template<typename C> static no  test(...);
  public:
    static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
};

//! \brief Detect if T defines function "expect::description"
//! \see http://stackoverflow.com/a/257382
template< typename T >
struct HasFunctionExpectDescription {
  private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;
    template<typename C> static yes test(decltype(&C::expect::description));
    template<typename C> static no  test(...);
  public:
    static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
};

//! \brief Detect if T defines function "expect::choices"
//! \see http://stackoverflow.com/a/257382
template< typename T >
struct HasFunctionExpectChoices {
  private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;
    template<typename C> static yes test(decltype(&C::expect::choices));
    template<typename C> static no  test(...);
  public:
    static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
};

} // tk::

#endif // Has_h
