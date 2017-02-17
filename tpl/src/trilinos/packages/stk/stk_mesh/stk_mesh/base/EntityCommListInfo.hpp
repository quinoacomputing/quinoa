#ifndef STK_ENTITYCOMMLIST_INFO_HPP
#define STK_ENTITYCOMMLIST_INFO_HPP

#include <vector>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk {
namespace mesh {

class Bucket;
struct EntityComm;

struct EntityCommListInfo
{
  EntityKey key;
  Entity    entity; // Might be invalid if entity has been deleted.
  Bucket* bucket;
  size_t bucket_ordinal;
  int  owner;
  const EntityComm* entity_comm; // Might be NULL if entity has been deleted.
};

inline
bool operator<(const EntityKey& key, const EntityCommListInfo& comm) { return key < comm.key; }

inline
bool operator<(const EntityCommListInfo& comm, const EntityKey& key) { return comm.key < key; }

inline
bool operator<(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs) { return lhs.key < rhs.key; }

inline
bool operator==(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs) { return lhs.key == rhs.key; }

inline
bool operator!=(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs) { return !(lhs == rhs); }

typedef std::vector<EntityCommListInfo> EntityCommListInfoVector;

struct IsInvalid
{
  bool operator()(const EntityCommListInfo& comm) const
  {
    return comm.key == EntityKey();
  }
};

}
}

#endif
