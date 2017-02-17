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

#ifndef __Teuchos_MatrixMarket_split_hpp
#define __Teuchos_MatrixMarket_split_hpp

#include <string>
#include <vector>


namespace Teuchos {
  namespace MatrixMarket {
    namespace details {

      //! Trim whitespace from both sides of the given string.
      std::string
      trim (const std::string& in);

      //! Return lowercase version of the given string.
      std::string
      lowercase (const std::string& in);

      //! Trim whitespace from both sides, and make lowercase.
      std::string
      trim_and_lowercase (const std::string& in);

      /// \brief Split the given string using the given set of delimiters.
      ///
      /// Split the string \c str, optionally starting at position \c
      /// start, into zero or more tokens separated by one or more of the
      /// given delimiter characters in \c delimiters.
      ///
      /// \param str [in] String to split into tokens
      /// \param delimiters [in] Array of one or more delimiter character(s)
      /// \param start [in] Position in \c str where the search should begin.
      ///   Defaults to zero.
      ///
      /// \return Vector of zero or more tokens, none of which contain any
      ///   of the delimiter character(s)
      std::vector<std::string>
      split (const std::string& str,
             const std::string& delimiters,
             const size_t start=0);

    } // namespace details
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_split_hpp
