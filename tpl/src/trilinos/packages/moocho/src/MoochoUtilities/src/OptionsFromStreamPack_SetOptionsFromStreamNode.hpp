// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef SET_OPTIONS_FROM_STREAM_NODE_H
#define SET_OPTIONS_FROM_STREAM_NODE_H

#include "OptionsFromStreamPack_SetOptionsFromStream.hpp"
#include "OptionsFromStreamPack_StringToIntMap.hpp"

namespace OptionsFromStreamPack {

/** \brief Node class for setting options from a stream.
  *
  * This class uses the template method pattern to
  * delegate the setting of options.
  */
class SetOptionsFromStreamNode: public SetOptionsFromStream {
public:

  /** \brief Constructs with the name of the options group and the names
    * of the options.
    *
    *	@param	options_group	The name of the options group to access
    *	@param	num_options		The number of options in the opitons
    *							group.
    *	@param	option_name		An array (length num_options) containing
    *							the names of the options.
    *	@param	exists_optional	Specifies if the options group must exist.
    */
  SetOptionsFromStreamNode( const std::string& options_group
    , int num_options, const char* option_names[]
    , bool exists_optional = true );

  /** \brief Overridden from SetOptionsFromStream and calls setOption(...).
    *
    * The options group #options_group# is used.  If this options
    * group does not exist and #exists_optional# == false then
    * an #std::invalid_argument# exception will be thrown.
    */
  void set_options( const OptionsFromStream& options );

protected:

  /** \brief To be overridden by the subclass to set an option given
    * its integer position and the option value.
    *
    * The integer possition returned is the possition of the option
    * in option_names[option_num] that was passed to the constructor.
    */
  virtual void setOption( int option_num, const std::string& option_value ) = 0;

private:
  StringToIntMap	name_map_;
  bool			exists_optional_;

};	// end class SetOptionsFromStreamNode

}	// end namespace OptionsFromStreamPack

#endif	// SET_OPTIONS_FROM_STREAM_NODE_H
