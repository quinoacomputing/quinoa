/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
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

#include "SundancePathUtils.hpp"
#include "SundanceDefaultPath.hpp"
#include <unistd.h>
#ifndef _MSC_VER
#include <sys/unistd.h>
#endif

using Teuchos::Array;

using std::ifstream;

namespace Sundance
{
string searchForFile(const std::string& name)
{
  std::string pathSep = "/";
  Array<string> path = parsePathStr();

  if (name.length() && name[0]=='/') return name; // Use absolute path!
  for (int i=0; i<path.size(); i++)
  {
    ifstream fileToTry((path[i] + pathSep + name).c_str());
    if (!fileToTry) continue;
    return path[i] + pathSep + name;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "could not find file "
    << name << " in path " << path);
}

string getPathStr() 
{
  char* pathEnvStr = getenv("SUNDANCE_PATH");
  char* pyPathEnvStr = getenv("PYTHONPATH");
  std::string path;
  
  if (pathEnvStr == NULL) 
  {
    path = defaultSundancePath();
  }
  else
  {
    path = pathEnvStr;
  }
  if (pyPathEnvStr!=NULL)
  {
    path = std::string(pyPathEnvStr) + ":" + path; 
  }
  return path;
}

Array<string> parsePathStr() 
{
  std::string pathStr = getPathStr();
  
  Array<string> rtn;

  unsigned int begin;
  unsigned int end;
  
  begin = pathStr.find_first_not_of(":");
  
  while (begin < pathStr.length())
  {
    end = pathStr.find_first_of(":", begin);

    rtn.append(pathStr.substr(begin, end-begin));
    begin = pathStr.find_first_not_of(":", end);
  }

  return rtn;
}
}
