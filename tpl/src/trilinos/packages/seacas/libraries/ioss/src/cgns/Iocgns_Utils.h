// Copyright(C) 2016
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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

#ifndef IOSS_IOCGNS_UTILS_H
#define IOSS_IOCGNS_UTILS_H

#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_Region.h>
#include <Ioss_SideBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_Utils.h>
#include <cgnslib.h>
#include <ostream>
#include <string>

namespace Iocgns {
  class Utils
  {
  public:
    Utils()  = default;
    ~Utils() = default;

    static void cgns_error(int cgnsid, const char *file, const char *function, int lineno,
                           int processor);

    template <typename INT>
    static void map_cgns_face_to_ioss(const Ioss::ElementTopology *parent_topo, size_t num_to_get,
                                      INT *idata)
    {
      // The {topo}_map[] arrays map from CGNS face# to IOSS face#.
      // See http://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid
      // NOTE: '0' for first entry is to account for 1-based face numbering.

      switch (parent_topo->shape()) {
      case Ioss::ElementShape::HEX:
        static int hex_map[] = {0, 5, 1, 2, 3, 4, 6};
        for (size_t i = 0; i < num_to_get; i++) {
          idata[2 * i + 1] = hex_map[idata[2 * i + 1]];
        }
        break;

      case Ioss::ElementShape::TET:
        static int tet_map[] = {0, 4, 1, 2, 3};
        for (size_t i = 0; i < num_to_get; i++) {
          idata[2 * i + 1] = tet_map[idata[2 * i + 1]];
        }
        break;

      case Ioss::ElementShape::PYRAMID:
        static int pyr_map[] = {0, 5, 1, 2, 3, 4};
        for (size_t i = 0; i < num_to_get; i++) {
          idata[2 * i + 1] = pyr_map[idata[2 * i + 1]];
        }
        break;

      case Ioss::ElementShape::WEDGE:
#if 0
	  static int wed_map[] = {0, 1, 2, 3, 4, 5}; // Same
	  // Not needed -- maps 1 to 1
	  for (size_t i=0; i < num_to_get; i++) {
	    idata[2*i+1] = wed_map[idata[2*i+1]];
	  }
#endif
        break;
      default:;
      }
    }

    static size_t resolve_nodes(Ioss::Region &region, int my_processor);

    static std::string map_cgns_to_topology_type(CG_ElementType_t type);
    static void add_sidesets(int cgnsFilePtr, Ioss::DatabaseIO *db);
    static void add_structured_boundary_conditions(int cgnsFilePtr, Ioss::StructuredBlock *block);
  };
}

#endif
