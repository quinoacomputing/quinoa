// *****************************************************************************
/*!
  \file      src/Control/Escaper.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     String escaper pulled over from PEGTL v0.32
  \details   String escaper pulled over from PEGTL v0.32 used to extract an
             std::string from pegtl::string.
*/
// *****************************************************************************
#ifndef Escaper_h
#define Escaper_h

#include <string>

namespace kw {

inline void escape_impl( std::string & result, const int i )
{
   switch ( i )
   {
      case '"':
         result += "\\\"";
         break;
      case '\\':
         result += "\\\\";
         break;
      case '\a':
         result += "\\a";
         break;
      case '\b':
         result += "\\b";
         break;
      case '\t':
         result += "\\t";
         break;
      case '\n':
         result += "\\n";
         break;
      case '\r':
         result += "\\r";
         break;
      case '\v':
         result += "\\v";
         break;
      case 32: case 33:          case 35: case 36: case 37: case 38: case 39:
      case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49:
      case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59:
      case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69:
      case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79:
      case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89:
      case 90: case 91:          case 93: case 94: case 95: case 96: case 97: case 98: case 99:
      case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109:
      case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119:
      case 120: case 121: case 122: case 123: case 124: case 125: case 126:
         result += char( i );
         break;
      default: {
         char tmp[ 12 ];
         ::snprintf( tmp, sizeof( tmp ), "\\u%04x", i );
         result += tmp;
      }  break;
   }
}

inline std::string escape( const int i )
{
   std::string nrv;
   escape_impl( nrv, i );
   return nrv;
}

inline std::string escape( const std::string & s )
{
   std::string nrv;

   for ( std::string::size_type i = 0; i < s.size(); ++i ) {
      escape_impl( nrv, s[ i ] );
   }
   return nrv;
}

template< int ... Cs > struct escaper;

template<>
struct escaper<>
{
   static std::string result()
   {
      return std::string();
   }
};

template< int C, int ... Cs >
struct escaper< C, Cs ... >
{
   static std::string result()
   {
      return escape( C ) + escaper< Cs ... >::result();
   }
};

} // ::kw

#endif // Escaper_h
