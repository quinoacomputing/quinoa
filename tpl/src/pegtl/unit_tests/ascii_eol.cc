// Copyright (c) 2016 Dr. Colin Hirsch and Daniel Frey
// Please see LICENSE for license or visit https://github.com/ColinH/PEGTL/

#include "test.hh"

namespace pegtl
{
   void unit_test()
   {
      verify_analyze< eol >( __LINE__, __FILE__, true, false );

      verify_rule< eol >( __LINE__, __FILE__,  "", result_type::LOCAL_FAILURE, 0 );

      for ( char i = 1; i < 127; ++i ) {
         verify_char< eol >( __LINE__, __FILE__, i, ( i == '\n' ) ? result_type::SUCCESS : result_type::LOCAL_FAILURE );
      }
      verify_rule< eol >( __LINE__, __FILE__,  "\r\n", result_type::SUCCESS, 0 );
      verify_rule< eol >( __LINE__, __FILE__,  "\n\r", result_type::SUCCESS, 1 );
      verify_rule< eol >( __LINE__, __FILE__,  "\na", result_type::SUCCESS, 1 );
      verify_rule< eol >( __LINE__, __FILE__,  "\ra", result_type::LOCAL_FAILURE, 2 );
      verify_rule< eol >( __LINE__, __FILE__,  "\r\na", result_type::SUCCESS, 1 );
      verify_rule< eol >( __LINE__, __FILE__,  "\r\n\r", result_type::SUCCESS, 1 );
      verify_rule< eol >( __LINE__, __FILE__,  "\r\n\n", result_type::SUCCESS, 1 );
   }

} // pegtl

#include "main.hh"
