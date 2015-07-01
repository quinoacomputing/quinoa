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


#include "SundanceFileIOChacoPartitioner.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;

using std::ofstream;
using std::ifstream;
using std::endl;


FileIOChacoPartitioner::FileIOChacoPartitioner(const std::string& filename,
  bool ignoreGhosts)
  : SerialPartitionerBase(ignoreGhosts), filename_(filename)
{}

void FileIOChacoPartitioner::writeGraph(const Mesh& mesh) const 
{
  Array<Array<int> > neighbors;
  int nEdges;

  getNeighbors(mesh, neighbors, nEdges);

  std::string gf = filename_ + ".graph";
  ofstream os(gf.c_str());

  os << neighbors.size() << " " << nEdges << std::endl;

  for (int i=0; i<neighbors.size(); i++)
  {
    for (int j=0; j<neighbors[i].size(); j++) 
    {
      if (j > 0) os << " ";
      os << neighbors[i][j]+1; // need unit offset here for Chaco
    }
    os << "\n";
  }
}


void FileIOChacoPartitioner::runChaco(int np) const 
{
  ofstream pf("User_Params");
  pf << 
    "OUTPUT_ASSIGN=true\n"
    "PROMPT=false\n"
    "ARCHITECTURE=1\n"
    "REFINE_PARTITION=4\n"
    "REFINE_MAP=true\n"
    "KL_BAD_MOVES=20\n"
    "KL_NTRIES_BAD=10\n"
    "KL_IMBALANCE=0.02\n"
    "INTERNAL_VERTICES=true\n"
    "MATCH_TYPE=2\n"
    "HEAVY_MATCH=true\n"
    "TERM_PROP=true\n"
    "COARSE_NLEVEL_KL=1\n"
    "COARSEN_RATIO_MIN=0.7\n"
    "CUT_TO_HOP_COST=1.0\n"
    "RANDOM_SEED=12345\n" << std::endl;

  ofstream chIn("chacoInput");
  chIn << filename_ + ".graph\n" << filename_ + ".assign\n1\n100\n"
       << np << "\n1\nn" << std::endl;

  int status = system("chaco < chacoInput");
  TEUCHOS_TEST_FOR_EXCEPTION(status < 0, std::runtime_error, "error detected in system call to run chaco");
}

bool FileIOChacoPartitioner::isEmptyLine(const std::string& x) const 
{
  return x.length()==0 || StrUtils::isWhite(x);
}

bool FileIOChacoPartitioner::getNextLine(std::istream& is, std::string& line,
                                         Array<string>& tokens,
                                         char comment) const 
{
  bool rtn = false;
  while ((rtn=StrUtils::readLine(is, line)))
    {
      if (line.length() > 0) line = StrUtils::before(line,comment);
      if (isEmptyLine(line)) continue;
      if (line.length() > 0) break;
    }
  tokens = StrUtils::stringTokenizer(line);
  return rtn;
}

void FileIOChacoPartitioner::getAssignments(const Mesh& mesh, int np,
  Array<int>& assignments) const 
{
  writeGraph(mesh);
  runChaco(np);

  std::string af = filename_ + ".assign";
  ifstream is(af.c_str());

  std::string line;
  Array<string> tokens;
    
  while (getNextLine(is, line, tokens, '#'))
  {
    assignments.append(StrUtils::atoi(tokens[0]));
  }
}

