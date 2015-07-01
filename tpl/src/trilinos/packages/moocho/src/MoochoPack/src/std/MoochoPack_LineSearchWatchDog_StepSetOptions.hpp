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

#ifndef LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
#define LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H

#include "MoochoPack_LineSearchWatchDog_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for LineSearchWatchDog_Step from a OptionsFromStream
  * object.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group LineSearchWatchDog {
    opt_kkt_err_threshold	= 1e-3; *** (+dbl)
    feas_kkt_err_threshold	= 1e-3; *** (+dbl)
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[opt_kkt_err_threshold] ToDo : Finish.
  *		Example: opt_kkt_err_threshold = 1e-1;
  *	\item[feas_kkt_err_threshold] ToDo : Finish.
  *		Example: feas_kkt_err_threshold = 1e-2;
  *	\end{description}
  */
class LineSearchWatchDog_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      LineSearchWatchDog_Step >
{
public:

  /** \brief . */
  LineSearchWatchDog_StepSetOptions(
    LineSearchWatchDog_Step* target = 0 );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class LineSearchWatchDog_StepSetOptions

}	// end namespace MoochoPack

#endif	// LINE_SEARCH_WATCH_DOG_STEP_SET_OPTIONS_H
