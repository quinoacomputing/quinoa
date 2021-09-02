// *****************************************************************************
/*!
  \file      src/Base/TeeBuf.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Tee stream buffer
  \details   Tee stream buffer that can be used to duplicate a stream. Example:
    \code{.cpp}
      // createa a file
      std::ofstream file( "data.txt" );
      // tbuf will write to both file and cout
      teebuf tbuf( file.rdbuf(), std::cout.rdbuf() );
      // replace cout's streambuf with tbuf
      scoped_streambuf_assignment ssa( std::cout, &tbuf );
      // write to both std::cout and data.txt
      std::cout << "Hello World" << std::endl;
   \endcode
*/
// *****************************************************************************
#ifndef Tee_h
#define Tee_h

#include <iostream>
#include <streambuf>
#include <fstream>
#include <memory>

namespace tk {

template< class charT, class traits = std::char_traits<charT> >
class basic_teebuf : public std::basic_streambuf< charT, traits > {

  public:
    using char_type = charT;
    using int_type = typename traits::int_type;
    using pos_type = typename traits::pos_type;
    using off_type = typename traits::off_type;
    using traits_type = traits;
    using streambuf_type = std::basic_streambuf< charT, traits >;

  private:
    streambuf_type* m_sbuf1;
    streambuf_type* m_sbuf2;
    std::unique_ptr< char_type[] > m_buffer;

    enum { BUFFER_SIZE = 4096 / sizeof( char_type ) };

  public:
    basic_teebuf( streambuf_type *sbuf1, streambuf_type *sbuf2 )
      : m_sbuf1( sbuf1 ), m_sbuf2( sbuf2 ),
        m_buffer( std::make_unique< char_type[] >( BUFFER_SIZE ) )
    { this->setp( m_buffer.get(), m_buffer.get() + BUFFER_SIZE ); }

    ~basic_teebuf() override { this->pubsync(); }

  protected:
    virtual int_type overflow( int_type c = traits_type::eof() ) override {
      // empty our buffer into m_sbuf1 and m_sbuf2
      std::streamsize n =
        static_cast< std::streamsize >( this->pptr() - this->pbase() );
      std::streamsize size1 = m_sbuf1->sputn( this->pbase(), n );
      std::streamsize size2 = m_sbuf2->sputn( this->pbase(), n );
      if ( size1 != n || size2 != n ) return traits_type::eof();

      // reset our buffer
      this->setp( m_buffer.get(), m_buffer.get() + BUFFER_SIZE );

      // write the passed character if necessary
      if ( !traits_type::eq_int_type(c, traits_type::eof()) ) {
        traits_type::assign( *this->pptr(), traits_type::to_char_type(c) );
        this->pbump(1);
      }

      return traits_type::not_eof(c);
    }

    virtual int sync() override {
      // flush our buffer into m_sbuf1 and m_sbuf2
      int_type c = this->overflow(traits_type::eof());

      // checking return for eof.
      if (traits_type::eq_int_type(c, traits_type::eof()))
          return -1;

      // flush m_sbuf1 and m_sbuf2
      if (m_sbuf1->pubsync() == -1 || m_sbuf2->pubsync() == -1)
          return -1;

      return 0;
    }
};

using teebuf = basic_teebuf< char >;
using wteebuf = basic_teebuf< wchar_t >;

template< class charT, class traits = std::char_traits< charT > >
struct scoped_basic_streambuf_assignment {
  using stream_type = std::basic_ios< charT, traits >;
  using streambuf_type = std::basic_streambuf< charT, traits >;

  stream_type& m_s;
  streambuf_type* m_orig_sb;

  scoped_basic_streambuf_assignment( stream_type &s, streambuf_type *new_sb )
    : m_s(s) { m_orig_sb = m_s.rdbuf( new_sb ); }

  ~scoped_basic_streambuf_assignment() { m_s.rdbuf( m_orig_sb ); }
};

using scoped_streambuf_assignment = scoped_basic_streambuf_assignment<char>;
using scoped_wstreambuf_assignment = scoped_basic_streambuf_assignment<wchar_t>;

} // tk::

#endif // TeeBuf_h
