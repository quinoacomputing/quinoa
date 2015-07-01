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


#include "SundanceMeshReaderBase.hpp"
#include "PlayaExceptions.hpp"
#include "SundancePathUtils.hpp"
#include "SundanceOut.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


MeshReaderBase::MeshReaderBase(const ParameterList& params)
  : MeshSourceBase(params), 
    filename_()
{
  filename_ = params.get<string>("Filename");
}



int MeshReaderBase::atoi(const std::string& x) const 
{
#ifndef TFLOP
  return std::atoi(x.c_str());
#else
  return ::atoi(x.c_str());
#endif
}

double MeshReaderBase::atof(const std::string& x) const 
{
#ifndef TFLOP
  return std::atof(x.c_str());
#else
  return ::atof(x.c_str());
#endif
}

bool MeshReaderBase::isEmptyLine(const std::string& x) const 
{
  return x.length()==0 || StrUtils::isWhite(x);
}

bool MeshReaderBase::getNextLine(std::istream& is, std::string& line,
                                         Array<string>& tokens,
                                         char comment) const 
{
  bool rtn = false;
  while ((rtn=StrUtils::readLine(is, line)))
    {
      SUNDANCE_OUT(this->verb() >= 3,
                   "read line [" << line << "]");

      if (line.length() > 0) line = StrUtils::before(line,comment);
      if (isEmptyLine(line)) continue;
      if (line.length() > 0) break;
    }
  tokens = StrUtils::stringTokenizer(line);
  return rtn;
}

RCP<std::ifstream> MeshReaderBase::openFile(const std::string& fname, 
                                               const std::string& description) const
{
  std::string f = searchForFile(fname);
  RCP<std::ifstream> rtn = rcp(new std::ifstream(f.c_str()));

  SUNDANCE_OUT(this->verb() > 2,
               "trying to open " << description << " file " << f);

  TEUCHOS_TEST_FOR_EXCEPTION(rtn.get()==0 || *rtn==0, std::runtime_error, 
                     "MeshReaderBase::openFile() unable to open "
                     << description << " file " << f);

  SUNDANCE_OUT(this->verb() > 0,
               "reading " << description << " from " << fname);

  return rtn;
}
