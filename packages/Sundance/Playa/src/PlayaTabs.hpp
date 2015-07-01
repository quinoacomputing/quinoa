/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef PLAYA_TABS_H
#define PLAYA_TABS_H

#include "PlayaDefs.hpp"

namespace Playa
{
/**
 * Tabbing utility for output. Constructing a new Tabs object automatically
 * increments the number of tabs to be written. When the Tabs object goes out
 * of scope, the original tabs level is restored. 
 *
 * The tab size and character can be specified through the setTabSize() and
 * setTabChar() methods, for example,
 * \code
 * Tabs::setTabChar('*');
 * Tabs::setTabSize(4);
 * \endcode
 * The tab character can be set on an object-by-object basis
 * through a constructor argument.
 *
 * By default, a header giving the depth of tabs is written to each line; this
 * can simplify scanning by eye for when a given tab level is reached. 
 * This header can be turned off by calling
 * \code
 * Tabs::showDepth() = false;
 * \endcode
 * 
 * Example: the code
 * \code
 * void f()
 * {
 *   Tabs tab;
 *   cout << tab << "in f()" << std::endl;
 *   g();
 *   cout << tab << "leaving f()" << std::endl;
 * }
 *
 * void g()
 * {
 *   Tabs tab0;
 *   cout << tab0 << "in g()" << std::endl;
 *   for (int i=0; i<3; i++)
 *     {
 *       Tabs tab1();
 *       cout << tab1 << "i=" << i << std::endl;
 *     }
 *   cout << tab0 << "leaving g()" << std::endl;
 * }
 * \endcode
 * writes the following output 
 * \code
 * [0]  in f()
 * [1]    in g()
 * [2]------i=0
 * [2]------i=1
 * [2]------i=2
 * [1]    leaving g()
 * [0]  leaving f()
 * \endcode 
 */
class Tabs
{
public:
  /** Constructor increments tab level */
  Tabs(bool jump=true);

  /** Destructor decrements tab level */
  ~Tabs();

  /** 
   * Print to stream. This method is usually not called directly, as
   * tabs will usually be written with the insertion operator
   */
  void print(std::ostream& os) const ;

  /** Change the tab size. Default is 2.  */
  static void setTabSize(int ts) {tabSize() = ts;}

  /** Indicate whether to print the tab depth as a header for each line. */
  static bool& showDepth() {static bool rtn = true; return rtn;}

private:
  /** */
  static int& tabLevel() {static int rtn = 0; return rtn;}

  /** */
  static int& tabSize() {static int rtn = 2; return rtn;}

  bool jump_;

  int myLevel_;
};
}

namespace Playa
{
/** \relates Tabs stream insertion operator for tab */
inline std::ostream& operator<<(std::ostream& os, const Tabs& t)
{
  t.print(os);
  return os;
}
}

#endif
