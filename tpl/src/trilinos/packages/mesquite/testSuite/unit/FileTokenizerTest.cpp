#include "Mesquite_FileTokenizer.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_cppunit/extensions/HelperMacros.h"
#include <assert.h>

#include <iostream>
using std::cout;

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef WIN32
#include <string>
#include <direct.h>
#endif

using Mesquite::FileTokenizer;
using Mesquite::MsqPrintError;
  
extern const char* const tokens[] = 
  { "foo", "bar", "0", "123abc", "123", "0x5", "1e200", "jason", "abc123", "1.0", 0 };
extern const char* const spaces[] = 
  { " ", " \n ", "\n\n\n ", "\t", "\t\t", "     ", "\n", "  ", "\t\n", 0 };

extern const bool doubles[]  = { false, false, true, false,  true, false,  true, false, false,  true };
extern const bool longs[]    = { false, false, true, false,  true,  true, false, false, false, false };
extern const bool booleans[] = { false, false, true, false, false, false, false, false, false, false };

class FileTokenizerTest : public CppUnit::TestFixture
{
private:

  CPPUNIT_TEST_SUITE( FileTokenizerTest );
  CPPUNIT_TEST( token_test );
  CPPUNIT_TEST( line_number_test );
  CPPUNIT_TEST( newline_test );
  CPPUNIT_TEST( match_one_test );
  CPPUNIT_TEST( match_multiple_test );
  CPPUNIT_TEST( double_test );
  CPPUNIT_TEST( long_test );
  CPPUNIT_TEST( boolean_test );
  CPPUNIT_TEST( unget_test );
  CPPUNIT_TEST_SUITE_END();
  
public:

  FileTokenizerTest() {}

  void setUp() { }
  
  void tearDown() { }
  
  FILE* make_file()
  {
#ifdef WIN32
    char name1[L_tmpnam_s];
    errno_t err = tmpnam_s( name1, L_tmpnam_s );
    
    // Get the current working directory: 
    char* buffer;
    buffer = _getcwd( NULL, 0 );
    std::string full_path(buffer);
    std::string temp_name(name1);
    full_path = full_path + temp_name;
    FILE* file = fopen( full_path.c_str(), "w+");
#else
    FILE* file = tmpfile();
#endif

    CPPUNIT_ASSERT( !!file );
    const char* const* t_iter = tokens;
    const char* const* s_iter = spaces;
    fputs( *t_iter, file ); ++t_iter;
    while (*t_iter)
    {
      fputs( *s_iter, file ); ++s_iter;
      fputs( *t_iter, file ); ++t_iter;
    }
    rewind( file );
    return file;
  }
    
  
  void token_test()
  {
    MsqPrintError err(cout);
    const char* token;
    
    FileTokenizer ft( make_file() );
    
    for (const char* const* t_iter = tokens; *t_iter; ++t_iter)
    {
      token = ft.get_string( err );
      CPPUNIT_ASSERT( token );
      CPPUNIT_ASSERT( !err );
      CPPUNIT_ASSERT( !strcmp( *t_iter, token ) );
    }
    
    token = ft.get_string(err);
    CPPUNIT_ASSERT( !token );
    CPPUNIT_ASSERT( ft.eof() );
    err.clear();
  }
  
  static int count_newlines( const char* str )
  {
    int result = 0;
    for (; *str; ++str)
      if (*str == '\n')
        ++result;
    return result;
  }
  
  void line_number_test()
  {
    MsqPrintError err(cout);
    const char* token;
    
    FileTokenizer ft( make_file() );
    
    int lineno = 1;
    token = ft.get_string( err );
    CPPUNIT_ASSERT( token );
    CPPUNIT_ASSERT( !err );
    
    for (const char* const* s_iter = spaces; *s_iter; ++s_iter)
    {
      token = ft.get_string( err );
      CPPUNIT_ASSERT( token );
      CPPUNIT_ASSERT( !err );
      
      lineno += count_newlines( *s_iter );
      CPPUNIT_ASSERT( ft.line_number() == lineno );
    }
  }
  
  void newline_test()
  {
    MsqPrintError err(cout);
    const char* token;
    
    FileTokenizer ft( make_file() );
    
    token = ft.get_string( err );
    CPPUNIT_ASSERT( token );
    CPPUNIT_ASSERT( !err );
    
    
    for (const char* const* s_iter = spaces; *s_iter; ++s_iter)
    {
      int count = count_newlines( *s_iter );
      bool b;
      while( count-- )
      {
        b = ft.get_newline( err );
        CPPUNIT_ASSERT( b && !err );
      }
      
      b = ft.get_newline( err );
      CPPUNIT_ASSERT( !b && err );
      err.clear();
    
      token = ft.get_string( err );
      CPPUNIT_ASSERT( token );
      CPPUNIT_ASSERT( !err );
    }
  }
  
  void match_one_test()
  {
    MsqPrintError err(cout);
    const char* const* t_iter;
    bool b;
    
    FileTokenizer ft( make_file() );
    
    for (t_iter = tokens; *t_iter; ++t_iter)
    {
      b = ft.match_token( *t_iter, err );
      CPPUNIT_ASSERT( b && !err );
    }
    
    FileTokenizer ft2( make_file() );
    
    b = ft2.match_token( "", err );
    CPPUNIT_ASSERT( !b && err );
    err.clear();
    
    b = ft2.match_token( "Mesquite", err );
    CPPUNIT_ASSERT( !b && err );
    err.clear();
  }
  
  void match_multiple_test()
  {
    MsqPrintError err(cout);
    const char* const* t_iter = tokens;
    
    FileTokenizer ft( make_file() );
    
    int result;
    const char* const test1[] = { *t_iter, "Mesquite", "x", 0 };
    result = ft.match_token( test1, err );
    CPPUNIT_ASSERT( result == 1 && !err );
    ++t_iter;
    
    const char* const test2[] = { "x", "y", *t_iter, 0 };
    result = ft.match_token( test2, err );
    CPPUNIT_ASSERT( result == 3 && !err );
    ++t_iter;
    
    const char* const test3[] = { *t_iter, 0 };
    result = ft.match_token( test3, err );
    CPPUNIT_ASSERT( result == 1 && !err );
    ++t_iter;
    
    const char* const test4[] = { "Mesquite", "Mesh", 0 };
    result = ft.match_token( test4, err );
    CPPUNIT_ASSERT( result == 0 && err );
    err.clear();
    ++t_iter;
  }
  
  void double_test()
  {
    MsqPrintError err(cout);
    FileTokenizer ft( make_file() );
    
    for (int i = 0; tokens[i]; ++i)
    {
      double value;
      ft.get_doubles( 1, &value, err );
      if (doubles[i])
      {
        CPPUNIT_ASSERT( !err );
        CPPUNIT_ASSERT( value == strtod( tokens[i], 0 ) );
      }
      else
      {
        CPPUNIT_ASSERT( err );
        err.clear();
      }
    }
  }
  
  void long_test()
  {
    MsqPrintError err(cout);
    FileTokenizer ft( make_file() );
    
    for (int i = 0; tokens[i]; ++i)
    {
      long value;
      ft.get_long_ints( 1, &value, err );
      if (longs[i])
      {
        CPPUNIT_ASSERT( !err );
        CPPUNIT_ASSERT( value == strtol( tokens[i], 0, 0 ) );
      }
      else
      {
        CPPUNIT_ASSERT( err );
        err.clear();
      }
    }
  }
   
  void boolean_test()
  {
    MsqPrintError err(cout);
    FileTokenizer ft( make_file() );
    
    for (int i = 0; tokens[i]; ++i)
    {
      bool value;
      ft.get_booleans( 1, &value, err );
      if (booleans[i])
      {
        CPPUNIT_ASSERT( !err );
        CPPUNIT_ASSERT( value == !!atoi( tokens[i] ) );
      }
      else
      {
        CPPUNIT_ASSERT( err );
        err.clear();
      }
    }
  }
   
  void unget_test()
  {
    MsqPrintError err(cout);
    FileTokenizer ft( make_file() );
    
    const char* const* t_iter = tokens;
    const char* token = ft.get_string(err);
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT( !strcmp( token, *t_iter ) );
    
    ft.unget_token();
    token = ft.get_string(err);
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT( !strcmp( token, *t_iter ) );
    
    ++t_iter;
    token = ft.get_string(err);
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT( !strcmp( token, *t_iter ) );
    
    ++t_iter;
    token = ft.get_string(err);
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT( !strcmp( token, *t_iter ) );
    
    ft.unget_token();
    token = ft.get_string(err);
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT( !strcmp( token, *t_iter ) );
  }    
   
};

    
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(FileTokenizerTest, "FileTokenizerTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(FileTokenizerTest, "Unit");


  
    
     
