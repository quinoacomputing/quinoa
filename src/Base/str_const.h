// *****************************************************************************
/*!
  \file      src/Base/str_const.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Constexpr, i.e., compile-time, string class
  \details   Constexpr, i.e., compile-time, string class
*/
// *****************************************************************************
#ifndef str_const_h
#define str_const_h

namespace tk {

//! \brief constexpr string
//! \author Scott Schurr
//! \see http://en.cppreference.com/w/cpp/language/constexpr
//! \see https://github.com/boostcon/cppnow_presentations_2012/blob/master/wed/schurr_cpp11_tools_for_class_authors.pdf?raw=true
class str_const {
  private:
    const char* const p_;
    const std::size_t sz_;
  public:
    // constructor
    template< std::size_t N >
    constexpr str_const( const char (&a)[N] ) : p_(a), sz_(N-1) {}
    // operator []
    constexpr char operator[] ( std::size_t n )
    { return n < sz_ ? p_[n] : throw std::out_of_range(""); }
    // size()
    constexpr std::size_t size() { return sz_; }
    // c_str accessor
    constexpr const char* c_str() const { return p_; }
};

} // tk::

#endif // str_const_h
