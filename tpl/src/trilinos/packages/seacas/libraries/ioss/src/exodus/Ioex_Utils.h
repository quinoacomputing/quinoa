#ifndef IOEX_UTILS_H
#define IOEX_UTILS_H
#include <Ioss_CoordinateFrame.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_Utils.h>

#include <assert.h>
#include <exodusII.h>
#include <set>
#include <string>
#include <vector>

// Contains code that is common between the file-per-processor (Iofx)
// and parallel exodus (Iopx) and base exodus (Ioex) classes.

namespace Ioss {
  class GroupingEntity;
  typedef std::vector<CoordinateFrame> CoordinateFrameContainer;
}

namespace Ioex {
  using EntityIdSet = std::set<std::pair<int64_t, int64_t>>;
  using SideSetSet  = std::set<std::string>;
  using SideSetMap  = std::map<std::string, const std::string, std::less<const std::string>>;

  struct TopologyMapCompare
  {
    bool operator()(const std::pair<std::string, const Ioss::ElementTopology *> &lhs,
                    const std::pair<std::string, const Ioss::ElementTopology *> &rhs) const
    {
      assert(lhs.second != nullptr);
      assert(rhs.second != nullptr);
      return lhs.first < rhs.first ||
             (!(rhs.first < lhs.first) && lhs.second->name() < rhs.second->name());
    }
  };

  using TopologyMap =
      std::map<std::pair<std::string, const Ioss::ElementTopology *>, int, TopologyMapCompare>;

  const char *Version();
  bool check_processor_info(int exodusFilePtr, int processor_count, int processor_id);

  void update_last_time_attribute(int exodusFilePtr, double value);
  bool read_last_time_attribute(int exodusFilePtr, double *value);

  bool type_match(const std::string &type, const char *substring);
  int64_t extract_id(const std::string &name_id);
  bool set_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset);
  int64_t get_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset);
  void decode_surface_name(Ioex::SideSetMap &fs_map, Ioex::SideSetSet &fs_set,
                           const std::string &name);
  void fix_bad_name(char *name);

  void exodus_error(int exoid, int lineno, const char *function, const char *filename);

  void check_non_null(void *ptr, const char *type, const std::string &name);

  int add_map_fields(int exoid, Ioss::ElementBlock *block, int64_t my_element_count,
                     size_t name_length);

  void add_coordinate_frames(int exoid, Ioss::Region *region);
  void write_coordinate_frames(int exoid, const Ioss::CoordinateFrameContainer &frames);

  inline int exodus_byte_size_api(int exoid)
  {
    // Check byte-size of integers stored on the database...
    int mode = ex_int64_status(exoid) & EX_ALL_INT64_API;
    if (mode) {
      return 8;
    }
    else {
      return 4;
    }
  }

  template <typename T> bool check_block_order(const std::vector<T *> &blocks)
  {
#ifndef NDEBUG
    // Verify that element blocks are defined in sorted offset order...
    typename std::vector<T *>::const_iterator I;

    int64_t eb_offset = -1;
    for (I = blocks.begin(); I != blocks.end(); ++I) {
      int64_t this_off = (*I)->get_offset();
      if (this_off < eb_offset) {
        {
          {
            return false;
          }
        }
      }
      eb_offset = this_off;
    }
#endif
    return true;
  }

  bool find_displacement_field(Ioss::NameList &fields, const Ioss::GroupingEntity *block, int ndim,
                               std::string *disp_name);

  char **get_exodus_names(size_t count, int size);
  void delete_exodus_names(char **names, int count);

  void get_fields(int64_t entity_count, char **names, size_t num_names,
                  Ioss::Field::RoleType fld_role, const char suffix_separator, int *local_truth,
                  std::vector<Ioss::Field> &fields);

  std::string get_entity_name(int exoid, ex_entity_type type, int64_t id,
                              const std::string &basename, int length, bool &db_has_name);

  void filter_element_list(Ioss::Region *region, Ioss::Int64Vector &elements,
                           Ioss::Int64Vector &sides, bool remove_omitted_elements);

  bool filter_node_list(Ioss::Int64Vector &               nodes,
                        const std::vector<unsigned char> &node_connectivity_status);

  template <typename T>
  void filter_node_list(T *data, std::vector<T> &dbvals,
                        const std::vector<int64_t> &active_node_index)
  {
    for (size_t i = 0; i < active_node_index.size(); i++) {
      data[i] = dbvals[active_node_index[i]];
    }
  }

  void filter_element_list(Ioss::Region *region, Ioss::Int64Vector &elements,
                           Ioss::Int64Vector &sides, bool remove_omitted_elements);

  void separate_surface_element_sides(Ioss::Int64Vector &element, Ioss::Int64Vector &sides,
                                      Ioss::Region *region, Ioex::TopologyMap &topo_map,
                                      Ioex::TopologyMap &    side_map,
                                      Ioss::SurfaceSplitType split_type,
                                      const std::string &    surface_name);
}
#endif
