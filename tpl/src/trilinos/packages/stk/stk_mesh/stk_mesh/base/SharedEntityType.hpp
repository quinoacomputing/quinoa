#ifndef STK_SHAREDENTITYTYPE_HPP
#define STK_SHAREDENTITYTYPE_HPP

#include <vector>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_topology/topology.hpp>

namespace stk {
namespace mesh {

struct shared_entity_type
{
  stk::topology::topology_t topology = stk::topology::INVALID_TOPOLOGY;
  std::vector<EntityKey>    nodes;
  stk::mesh::EntityKey      local_key;
  stk::mesh::EntityKey      global_key;
  std::vector<int>          sharing_procs;
  stk::mesh::Entity         entity;
  bool                      need_update_nodes = false;

  shared_entity_type(stk::mesh::EntityKey _key, stk::mesh::Entity _entity, stk::topology::topology_t _top) :
      topology(_top), local_key(_key), global_key(_key), entity(_entity) {}

  friend inline bool operator < (shared_entity_type const& l, shared_entity_type const& r)
  {
    if (l.topology < r.topology)   return true;
    if (l.topology > r.topology)   return false;
    return l.nodes < r.nodes;
  }

  friend inline bool operator == (shared_entity_type const& l, shared_entity_type const& r)
  {

    bool sameTopologyAndNodes =  (l.topology == r.topology) && (l.nodes==r.nodes);
    return sameTopologyAndNodes;
  }
};

}
}

#endif
