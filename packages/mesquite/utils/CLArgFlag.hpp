/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file CLArgFlag.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_CLARG_FLAG_HPP
#define MSQ_CLARG_FLAG_HPP

#include "Mesquite.hpp"
#include "CLArgs.hpp"
#include <string>
#include <vector>
#include <iosfwd>

class CLArgFlag 
{
  private:
    char mFlag;
    std::string mDesc;
  
  protected:  
    CLArgFlag( char flag, const char* desc )
      : mFlag(flag), mDesc(desc)
      {}

  public:
    virtual ~CLArgFlag();
  
    char flag() const { return mFlag; }
    const char* desc() const { return mDesc.c_str(); }
    virtual const CLArgs::ArgIBase* callback() const = 0;
  
    virtual bool is_toggle() const;
    virtual bool parse( const char* string_value ) const = 0;
    
      /** Get brief format of option */
    virtual std::string brief() const = 0;
      /** Get UNIX man-page formated synposis of flag */
    virtual std::string manstr() const = 0;

    virtual bool add_set( int size, const char* const* names );
    
    std::string make_man_string( int count, const char* names[] ) const;
    std::string make_man_string( const char* arg_name ) const
      { return make_man_string( 1, &arg_name ); }
    std::string make_man_string( ) const
      { return make_man_string( 0, 0 ); }
    std::string make_literal_man_string( int count, const char* literal_args[] ) const;
    std::string make_literal_man_string( const char* literal_args ) const
      { return make_literal_man_string( 1, &literal_args ); }
};

class CLArgToggle : public CLArgFlag
{
  private:
    bool mValue; //!< value to pass to callback when flag is encountered
    CLArgs::ToggleArgI* mCallback;
    CLArgToggle* mOpposite;
  public:
    CLArgToggle( char flag, 
                 const char* desc,
                 bool value,
                 CLArgs::ToggleArgI* callback )
      : CLArgFlag( flag, desc ),
        mValue( value ),
        mCallback( callback ),
        mOpposite( 0 )
      { }
    CLArgToggle( char flag, 
                 const char* desc,
                 CLArgToggle* opposite )
      : CLArgFlag( flag, desc ),
        mValue( !opposite->mValue ),
        mCallback( opposite->mCallback ),
        mOpposite( opposite )
      {
        mOpposite->mOpposite = this;
      }

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool is_toggle() const
      { return true; }
    virtual bool parse( const char* option ) const;
      
    virtual std::string brief() const;
    virtual std::string manstr() const;
};

class CLArgString : public CLArgFlag
{
  private:
    std::string mName;
    CLArgs::StringArgI* mCallback;
  
  public:
    CLArgString( char fl, const char* name, const char* desc,
                 CLArgs::StringArgI* callback )
                : CLArgFlag( fl, desc ),
                  mName( name ),
                  mCallback( callback )
                  {}

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool parse( const char* option ) const;
      

    virtual std::string brief() const;
    virtual std::string manstr() const;
};

class CLArgLong : public CLArgFlag
{
  private:
    CLArgs::LongArgI* mCallback;
    std::string mName;
  
  public:
    CLArgLong( char fl, const char* name, const char* desc,
               CLArgs::LongArgI* callback )
      : CLArgFlag( fl, desc ),
        mCallback( callback ),
        mName( name )
    {}

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool parse( const char* option ) const;

    virtual std::string brief() const;
    virtual std::string manstr() const;
};

class CLArgInt : public CLArgFlag
{
  private:
    CLArgs::IntArgI* mCallback;
    std::string mName;
  
  public:
    CLArgInt(  char fl, const char* name, const char* desc,
               CLArgs::IntArgI* callback )
      : CLArgFlag( fl, desc ),
        mCallback( callback ),
        mName( name )
    {}

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool parse( const char* option ) const;

    virtual std::string brief() const;
    virtual std::string manstr() const;
};

class CLArgDouble : public CLArgFlag
{
  private:
    CLArgs::DoubleArgI* mCallback;
    std::string mName;
  
  public:
    CLArgDouble( char fl, const char* name, const char* desc,
                 CLArgs::DoubleArgI* callback )
      : CLArgFlag( fl, desc ),
        mCallback( callback ),
        mName( name )
    {}

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool parse( const char* option ) const;

    virtual std::string brief() const;
    virtual std::string manstr() const;
};

  
class CLArgIDList : public CLArgFlag
{
  public:
    CLArgs::IntListArgI* mCallback;
    
  public:
    CLArgIDList( char fl, const char* desc, CLArgs::IntListArgI* callback )
      : CLArgFlag( fl, desc ), mCallback( callback )
      {}

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool parse( const char* str ) const;

    virtual std::string brief() const;
    virtual std::string manstr() const;
};

class CLArgListData 
{
  std::vector< std::vector< std::string > > mSets;

  public:
  
    bool add_set( int size, const char* const* names );
    bool acceptable_length( unsigned len ) const;
    bool accept_any_length() const
      { return mSets.empty(); }
    
    std::string set_string( int set ) const;
    std::string brief() const;
    std::string manstr( char type_char, const CLArgFlag& f ) const;
};

class CLArgIntList : public CLArgFlag
{
  private:
    CLArgListData listData;
    CLArgs::IntListArgI* mCallback;
  
  public:
    CLArgIntList( char fl, const char* desc, CLArgs::IntListArgI* callback )
      : CLArgFlag( fl, desc ), mCallback( callback )
    {}

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool parse( const char* option ) const;
      
    virtual std::string brief() const;
    virtual std::string manstr() const;
    virtual bool add_set( int size, const char* const* names )
      { return listData.add_set( size, names ); }
};

class CLArgDoubleList : public CLArgFlag
{
  private:
    CLArgListData listData;
    std::string mName;
    CLArgs::DoubleListArgI* mCallback;
  
  public:
    CLArgDoubleList( char fl, 
                     const char* desc,
                     CLArgs::DoubleListArgI* callback )
      : CLArgFlag( fl, desc ),
        mCallback( callback )
        {}

    virtual const CLArgs::ArgIBase* callback() const { return mCallback; }
    
    virtual bool parse( const char* option ) const;

    virtual std::string brief() const;
    virtual std::string manstr() const;
    virtual bool add_set( int size, const char* const* names )
      { return listData.add_set( size, names ); }
};

#endif
