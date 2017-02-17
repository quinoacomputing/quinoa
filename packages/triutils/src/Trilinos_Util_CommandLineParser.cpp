// @HEADER
// ***********************************************************************
//
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

/* ------------------------------------------------------------------------ */

#include "Trilinos_Util.h"
#include "Trilinos_Util_CommandLineParser.h"
#include <cstring>

using namespace std;

Trilinos_Util_Map::Trilinos_Util_Map(void)
{

  SetLabel("Trilinos_Util_Map");

  return;

}

void Trilinos_Util_Map::Reset(void)
{

  SetLabel("");
  /*FIXME
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    ci->pop_back();
  }
*/
  return;

}

// ================================================ ====== ==== ==== == =

Trilinos_Util::CommandLineParser::CommandLineParser(int argc, char *argv[])
{

  SetLabel("Trilinos_Util::CommandLineParser");

  char str[80];
  string value, param;

  Set("PROGRAM_NAME_",argv[0]);

  sprintf(str,"%d",argc);
  Set("_N_ARGS_",str);

  // first, manage possible arguments without specifier
  // (e.g., a.out 12 -nx 10 -ny 20)
  // Those arguments are called _ARGV_1_, ... , _ARGV_N_
  // and the value of N is given by _N_UNNAMED_ARGS_

  int N_args = 0;
  int i=1;

  for( i=1 ; i<argc ; ++i ) {
    if( *(argv[i]) == '-' ) break;
    N_args++;
    sprintf( str, "ARGV_%d", N_args );
    string param3;
    param3 = argv[i];
    Set(param3,value);
  }

  sprintf(str,"%d",N_args);
  Set("_N_UNNAMED_ARGS_",str);

  // now only arguments with a dash (possibly followed by one
  // other specifier)

  for( ; i<argc ; ++i ) {
    // check if the option has a `=' inside.
    // If so, split the string into two substrings
    char * pos = strchr( argv[i], '=');
    if( pos != NULL ) {
      *pos = '\0';
      param = argv[i], value = pos+1;
      Set(param,value);
    } else if( i<argc-1 ) {
      if( *(argv[i+1]) != '-' ) {
  param = argv[i], value = argv[i+1];
  Set(param,value);
  ++i;
      } else {
  param = argv[i], value = "";
  Set(param,value);
      }
    } else {
      param = argv[i], value = "";
      Set(param,value);
    }

  }

}

// ================================================ ====== ==== ==== == =

int Trilinos_Util_Map::Get( const string input, const int def_value)
{

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input )
      return( atoi(Map_[input].c_str()) );
  }

  return def_value;

}

// ================================================ ====== ==== ==== == =

double Trilinos_Util_Map::Get( const string input, const double def_value)
{

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input )
      return( atof(Map_[input].c_str()) );
  }

  return def_value;

}

// ================================================ ====== ==== ==== == =

string Trilinos_Util_Map::Get( const string input, const string def_value)
{

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input )
      return( Map_[input] );
  }

  return def_value;

}

// ================================================ ====== ==== ==== == =

bool Trilinos_Util_Map::Has( const string input)
{

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input )
      return true;
  }
  return false;

}

// ================================================ ====== ==== ==== == =

void Trilinos_Util_Map::ShowAll() const
{

  cout << "\n" << Label_ << " :: \n";

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first.at(0) != '_' )
      cout << (*ci).first << " = " << (*ci).second << endl;
  }

} /* ShowAll */

void Trilinos_Util_Map::ShowReallyAll() const
{

  cout << "\nTrilinos_Util_CommandLineParser :: \n";

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    cout << (*ci).first << " = " << (*ci).second << endl;
  }

} /* ShowReallyAll */

bool Trilinos_Util_Map::Add( const string input, const string value )
{

  // check that "input" has not been already inserted
  if( this->Has(input) == true )
    return false;

  Map_[input] = value;

  return true;

} /* AddOption */

bool Trilinos_Util_Map::Set( const string input, const double value )
{

  char value2[80];
  sprintf( value2, "%e", value);
  return( Set(input,value2) );
}

bool Trilinos_Util_Map::Set( const string input, const int value )
{

  char value2[80];
  sprintf( value2, "%d", value);
  return( Set(input,value2) );
}

bool Trilinos_Util_Map::Set( const string input, const string value )
{

  Map_[input] = value;

  return true;

} /* Set */

bool Trilinos_Util_Map::Set( const string input, const char * value )
{

  string val(value);

  Map_[input] = val;

  return true;

} /* Set */

string Trilinos_Util::CommandLineParser::GetProgramName( void )
{
  return( Get("_PROGRAM_NAME_", "UNDEFINED" ) );

}

int Trilinos_Util::CommandLineParser::GetIntShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0;

} /* GetIntShellVariable */

double Trilinos_Util::CommandLineParser::GetDoubleShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0.0;

} /* GetDoubleShellVariable */

string Trilinos_Util::CommandLineParser::GetStringShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer == NULL )
    return( "" );

  return( buffer );

} /* GetCharShellVariable */

// ================================================ ====== ==== ==== == =

ostream & operator << (ostream & os,
           const Trilinos_Util_Map & S)
{
  S.ShowAll();
  return os;
}

// ================================================ ====== ==== ==== == =

Trilinos_Util::InputFileReader::InputFileReader(const char FileName[]) :
  FileName_(FileName), CommentChars_("#"), SeparationChars_("="),
  FileHasBeenRead_(false)
{
}

Trilinos_Util::InputFileReader::~InputFileReader()
{

  FileName_ = "";
  CommentChars_ = "";
  SeparationChars_ = "";
  Reset();
  FileHasBeenRead_ = false;

}

string Trilinos_Util::InputFileReader::GetFileName() const
{
  return FileName_;
}

void Trilinos_Util::InputFileReader::SetCommentChars(const string c)
{
  CommentChars_ = c;
  return;
}

void Trilinos_Util::InputFileReader::SetSeparationChars(const string c)
{
  SeparationChars_ = c;
  return;
}

#include <iostream>
#include <fstream>

int Trilinos_Util::InputFileReader::ReadFile(const char * FileName)
{
  FileName_ = FileName;

  return( ReadFile() );
}

int Trilinos_Util::InputFileReader::ReadFile()
{

  ifstream File(FileName_.c_str());

  if( File.good() == false ) {
    std::cerr << "Error opening file `" << FileName_ << "'\n";
    return -1;
  }

  const int CharMax = 255;
  char line[CharMax];
  string Option, Value;

  while( File.eof() == false ) {

    File.getline(line,255);
    string StrLine = line;
    for( int k=0 ; k<(int)CommentChars_.length() ; ++k ) {
      int CommentPos = StrLine.find(CommentChars_.at(k));
      if( CommentPos != -1 ) {
  StrLine = StrLine.substr(0,CommentPos);
      }
    }
    int Length = StrLine.length();
    for( int k=0 ; k< (int) SeparationChars_.length() ; ++k ) {
      int SepPos = StrLine.find(SeparationChars_.at(k));
      if( SepPos > 0 ) {
  Option = StrLine.substr(0,SepPos);
  Value = StrLine.substr(SepPos+1,Length);
  // ~!@ to erase spaces...
  if( Option.length() > 0 ) Set(Option,Value);
  break;
      }
    }
  }

  // close file
  File.close();

  return 0;

}
