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

#ifndef OPTIONS_FORM_STEAM_PACK_COMMANDLINE_OPTIONS_FROM_STREAM_PROCESSOR_HPP
#define OPTIONS_FORM_STEAM_PACK_COMMANDLINE_OPTIONS_FROM_STREAM_PROCESSOR_HPP

#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace OptionsFromStreamPack {

/** \brief Reads from a file and/or parses from the commandline to initalize
 * an <tt>OptionsFromStream</tt> object.
 *
 * ToDo: Finish documentation!
 */
class CommandLineOptionsFromStreamProcessor {
public:

  /** \brief Construct with default values for the options file name and extra
   * options.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->get_options_file_name()==options_file_name</tt>
   * <li><tt>this->get_extra_options_str()==extra_options_str</tt>
   * </ul>
   *
   * <b>Suggestion:</b> Do not pass in arguments to the consturctor.  Instead,
   * construct the object and then call the setup frunction.
   */
  CommandLineOptionsFromStreamProcessor(
    const std::string  &options_file_name_opt_name   = "ofs-options-file"
    ,const std::string  &options_file_name_opt_doc   = "The name of the file containing input options for OptionsFromStream object."
    ,const std::string &options_file_name            = ""
    ,const std::string &extra_options_str_opt_name   = "ofs-extra-options"
    ,const std::string &extra_options_str_opt_doc    = "Extra options in format \"OptGrp1{name1=val1,...,namen=valn}:OptGr2{name1=val1,...,namen=valn}:...\""
    ,const std::string &extra_options_str            = ""
    );
  // RAB: 2006/01/27: Note, this value contains no semi-columns since this
  // conflicts with Trilinos' runtests script.  Therefore, I have to replace
  // the ',' separators with ';' below!
  // Note: we can leave off the last ',' since it turns out that the
  // way the new OptionsFromStream::parse_options(...) is written that
  // the last semicolon in an options group is not necessary!

  /** \brief Set the <tt>OptionsFromStream</tt> object that will be used to
   * fill options in to.
   *
   * Postconditions:<ul>
   * <li><tt>this->get_options()==options.get()</tt>
   * </ul>
   */
  void set_options(
    Teuchos::RCP<OptionsFromStream> const& options
    );

  /** \brief Just return the <tt>OptionsFromStream</tt> object in its current
   * state.
   *
   * This function does not force the options to be processed.
   */
  Teuchos::RCP<OptionsFromStream> get_options() const;

  /** \brief Set the name of the commandline option name that specifies the
   * options file name.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->options_file_name_opt_name()==options_file_name_opt_name</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,options_file_name_opt_name);

  /** \brief Set the documentation of the commandline option name that specifies the
   * options file name.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->options_file_name_opt_doc()==options_file_name_opt_doc</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,options_file_name_opt_doc);

  /** \brief Set the options file name manually (which will be used for
   * default value for "--${options_file_name_opt_name}=???" commandline option).
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->options_file_name()==options_file_name</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,options_file_name);

  /** \brief Set the name of the commandline option name that specifies the
   * extra options.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->extra_options_str_opt_name()==extra_options_str_opt_name</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,extra_options_str_opt_name);

  /** \brief Set the documentation of the commandline option name that
   * specifies the extra options.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->oextra_options_str_opt_doc()==extra_options_str_opt_doc</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,extra_options_str_opt_doc);

  /** \brief Set the extra commandline options string (can be used for default value
   * for "--${extra_options_str_opt_name}=???" commandline option).
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->get_extra_options_str()===extra_options_str</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS(std::string,extra_options_str);

  /** \brief Setup a comandline processor before it processes commandline
   * options or reads form a file.
   *
   * Sets up two options on the commandline processor:<ul>
   * <li><tt>--${this->options_file_name_opt_name()}</tt>: Sets the value returned from <tt>this->options_file_name()</tt>
   * <li><tt>--${this->extra_options_str_opt_name()}</tt>: Sets the value returned from <tt>this->extra_options_str()</tt>
   * </ul>
   */
  void setup_commandline_processor(
    Teuchos::CommandLineProcessor *clp
    );

  /** \brief Read the options file and/or process the commandline options.
   *
   * If <tt>this->options_file_name()!=""</tt>, then this file is read and
   * processed first.  Then, if <tt>this->extra_options_str()!=""</tt>, this
   * is followed by the processing of these extra options which may overide
   * options set in the input options file.
   *
   * <b>Note:</b> If <tt>this->get_options().get()==NULL</tt> before this
   * function is called and if <tt>this->options_file_name()!=""</tt> or
   * <tt>this->extra_options_str()!=""</tt>, the an <tt>OptionsFromStream</tt>
   * object will be created anew to store the parsed options.
   *
   * <b>Note:</b> Calling this function will result in the options file and
   * extra options to be parsed over and over again and add to or override the
   * options currently already in <tt>*this->get_options()</tt>.
   */
  void process_options();

  /** \brief Calls <tt>process_options()</tt> and returns
   * <tt>get_options()</tt>
   */
  Teuchos::RCP<OptionsFromStream> process_and_get_options();

private:

  Teuchos::RCP<OptionsFromStream>     options_;

};

}	// end namespace OptionsFromStreamPack

#endif	// OPTIONS_FORM_STEAM_PACK_COMMANDLINE_OPTIONS_FROM_STREAM_PROCESSOR_HPP
