/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_MATRIX_IO_DECL
#define TPETRA_MATRIX_IO_DECL

#include <string>
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Tpetra {
  namespace Utils {

    bool parseIfmt(ArrayRCP<char> fmt, int &perline, int &width);
    bool parseRfmt(ArrayRCP<char> fmt, int &perline, int &width, int &prec, char &flag);
    void readHBInfo(const std::string &filename, int &M, int &N, int &nz, ArrayRCP<char> &Type, int &Nrhs);

    void readHBHeader(std::ifstream &in_file, ArrayRCP<char> &Title, ArrayRCP<char> &Key, ArrayRCP<char> &Type, 
        int &Nrow, int &Ncol, int &Nnzero, int &Nrhs,
        ArrayRCP<char> &Ptrfmt, ArrayRCP<char> &Indfmt, ArrayRCP<char> &Valfmt, ArrayRCP<char> &Rhsfmt, 
        int &Ptrcrd, int &Indcrd, int &Valcrd, int &Rhscrd, ArrayRCP<char> &Rhstype);

    void readHBMatDouble(const std::string &filename, int &M, int &N, int &nonzeros, std::string &Type, ArrayRCP<int> &colptr, ArrayRCP<int> &rowind, ArrayRCP<double> &val);

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    void
    generateMatrix(const RCP<ParameterList> &plist,
                   const RCP<const Comm<int> > &comm, 
                   const RCP<Node> &node,
                   RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A);

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    void
    readHBMatrix(const std::string &filename, 
                 const RCP<const Comm<int> > &comm, 
                 const RCP<Node> &node,
                 RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A,
                 RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = null, 
                 const RCP<ParameterList> &params = null);

  } // end of Tpetra::Utils namespace
} // end of Tpetra namespace

#endif
