// *****************************************************************************
/*!
  \file      src/Base/Make_unique.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Define make_unique for unique_ptr until C++14
  \details   Define make_unique for unique_ptr until C++14. When we switch to
    C++14, this can go away. The code below is lifted from
    http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen, generated on
    2014-01-06.
*/
// *****************************************************************************
#ifndef Make_unique_h
#define Make_unique_h

#include <memory>

namespace tk {

// libc++ has make_unique in C++1y mode, but only for a commit later than
// r181765 (e.g., r203847 already has it), until that we'll use the one taken
// from gcc below

#if not defined (_LIBCPP_STD_VER) || (_LIBCPP_STD_VER <= 11)

using std::unique_ptr;
using std::remove_extent;

template<typename _Tp>
  struct _MakeUniq
  { typedef unique_ptr<_Tp> __single_object; };

template<typename _Tp>
  struct _MakeUniq<_Tp[]>
  { typedef unique_ptr<_Tp[]> __array; };

template<typename _Tp, size_t _Bound>
  struct _MakeUniq<_Tp[_Bound]>
  { struct __invalid_type { }; };

/// std::make_unique for single objects
template<typename _Tp, typename... _Args>
  inline typename _MakeUniq<_Tp>::__single_object
  make_unique(_Args&&... __args)
  { return unique_ptr<_Tp>(new _Tp(std::forward<_Args>(__args)...)); }

/// std::make_unique for arrays of unknown bound
template<typename _Tp>
  inline typename _MakeUniq<_Tp>::__array
  make_unique(size_t __num)
  { return unique_ptr<_Tp>(new typename remove_extent<_Tp>::type[__num]()); }

/// Disable std::make_unique for arrays of known bound
template<typename _Tp, typename... _Args>
  inline typename _MakeUniq<_Tp>::__invalid_type
  make_unique(_Args&&...) = delete;

#else

using std::make_unique;

#endif

} // tk::

#endif // Make_unique_h
