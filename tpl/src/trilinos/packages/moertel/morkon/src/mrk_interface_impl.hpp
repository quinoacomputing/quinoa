/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2015) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef MORKON_EXP_API_INTERFACE_IMPL_H
#define MORKON_EXP_API_INTERFACE_IMPL_H

namespace morkon_exp {

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
Interface<DeviceType, DIM,  FACE_TYPE >::Interface(Morkon_Manager<DeviceType, DIM, FACE_TYPE> * manager)
  : m_manager(manager)
  , m_committed(false)
{
  m_hs_adapters.resize(2, 0);
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
Interface<DeviceType, DIM,  FACE_TYPE >::~Interface()
{
  delete m_hs_adapters[0];
  delete m_hs_adapters[1];
  m_hs_adapters.clear();
}


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Interface<DeviceType, DIM,  FACE_TYPE >::hsa_add_node(SideEnum which_side, global_idx_t gbl_node_id, const double coords[])
{
  if (m_committed)
  {
    return false;
  }

  Interface_HostSideAdapter<DIM> *ifc_hsa = m_hs_adapters[which_side];
  if (!ifc_hsa)
  {
    ifc_hsa = m_hs_adapters[which_side] = new Interface_HostSideAdapter<DIM>();
  }

  typename Interface_HostSideAdapter<DIM>::node_map_type::iterator probe = ifc_hsa->m_nodes.find(gbl_node_id);
  if (probe != ifc_hsa->m_nodes.end())
  {
    return false;
  }

  typename Interface_HostSideAdapter<DIM>::NodeInfo node_info;
  node_info.m_id        = gbl_node_id;
  node_info.m_side      = which_side;
  node_info.m_coords[0] = coords[0];
  node_info.m_coords[1] = coords[1];
  node_info.m_coords[2] = (DIM == 2 ? 0 : coords[2]);
  ifc_hsa->m_nodes.insert(probe, std::make_pair(gbl_node_id, node_info));

  return true;
}


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Interface<DeviceType, DIM, FACE_TYPE>::hsa_add_face(SideEnum which_side, global_idx_t gbl_face_id, int num_nodes, const global_idx_t gbl_node_id[])
{
  if (m_committed )
  {
    return false;
  }

  Interface_HostSideAdapter<DIM> *ifc_hsa = m_hs_adapters[which_side];
  if (!ifc_hsa)
  {
    ifc_hsa = m_hs_adapters[which_side] = new Interface_HostSideAdapter<DIM>();
  }

  for (size_t node_i = 0; node_i < TopoConsts<FACE_TYPE>::NODES_PER_FACE; ++node_i)
  {
    typename Interface_HostSideAdapter<DIM>::node_map_type::iterator probe_node
      = ifc_hsa->m_nodes.find(gbl_node_id[node_i]);
    if (probe_node == ifc_hsa->m_nodes.end())
    {
      return false;
    }
  }

  typename Interface_HostSideAdapter<DIM>::face_map_type::iterator probe_face = ifc_hsa->m_faces.find(gbl_face_id);
  if (probe_face != ifc_hsa->m_faces.end())
  {
    return false;
  }

  typename Interface_HostSideAdapter<DIM>::FaceInfo face_info;
  face_info.m_id        = gbl_face_id;
  face_info.m_side      = which_side;
  face_info.m_nodes.resize(num_nodes);
  for (int i =  0; i < num_nodes; ++i)
  {
    face_info.m_nodes[i] = gbl_node_id[i];
  }

  ifc_hsa->m_faces.insert(probe_face, std::make_pair(gbl_face_id, face_info));

  return true;
}


}

#endif
