// Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef EXO_BLOCK_H
#define EXO_BLOCK_H

#include "exo_entity.h"
#include <iostream>

#include <string>

template <typename INT> class ExoII_Read;

template <typename INT> class Exo_Block : public Exo_Entity
{
public:
  Exo_Block();
  Exo_Block(int exo_file_id, size_t block_id);
  Exo_Block(int exo_file_id, size_t block_id, const char *type, size_t num_elmts,
            size_t num_nodes_per_elmt_x);
  ~Exo_Block() override;

  std::string Load_Connectivity();
  std::string Free_Connectivity();

  // Access functions:

  const std::string &Elmt_Type() const { return elmt_type; }
  size_t             Num_Nodes_per_Elmt() const { return num_nodes_per_elmt; }

  // Block description access functions:

  const INT *Connectivity() const { return conn; }  // 1-offset connectivity
  const INT *Connectivity(size_t elmt_index) const; // 1-offset connectivity

  std::string Give_Connectivity(size_t &num_e,      // Moves connectivity matrix
                                size_t &npe,        // to conn pointer and sets
                                INT *&  recv_conn); // its own to null.

  // Misc:

  int Check_State() const;

  void Display_Stats(std::ostream & = std::cout) const;
  void Display(std::ostream & = std::cout) const;

private:
  Exo_Block(const Exo_Block &);                  // Not written.
  const Exo_Block &operator=(const Exo_Block &); // Not written.

  void entity_load_params() override;

  EXOTYPE     exodus_type() const override;
  const char *label() const override { return "Element Block"; }
  const char *short_label() const override { return "block"; }

  std::string elmt_type;
  int         num_nodes_per_elmt;

  INT *conn; // Array; holds a matrix, num_elmts by num_nodes_per_elmt.

  friend class ExoII_Read<INT>;
};

#endif
