#ifndef READWRITESIDESETTESTER_HPP_
#define READWRITESIDESETTESTER_HPP_

#include <vector>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/InputFile.hpp>   // for InputFile for m_input_files
#include <stk_util/parallel/Parallel.hpp> // for stk::parallel_machine_size
#include "Ioss_Region.h"                // for Region, NodeSetContainer, etc
#include "Ioss_SideBlock.h"             // for SideBlock
#include "Ioss_SideSet.h"               // for SideSet, SideBlockContainer
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"

namespace stk{ namespace unit_test_util{ namespace sideset{

struct IdAndSideSet {
    int id;
    stk::mesh::SideSet sideSet;
};

typedef std::vector<IdAndSideSet> SideSetData;



enum ReadMode { READ_SERIAL_AND_DECOMPOSE, READ_ALREADY_DECOMPOSED };

class StkMeshIoBrokerTester : public stk::io::StkMeshIoBroker
{
public:
    StkMeshIoBrokerTester() {}

    virtual void populate_mesh(bool delay_field_data_allocation = true)
    {
        stk::io::StkMeshIoBroker::populate_mesh(delay_field_data_allocation);
        extract_sideset_data_from_io(sideset_data);
    }

private:
    int convert_to_zero_based_ordinal(int one_based_ordinal) const
    {
        return one_based_ordinal-1;
    }

    stk::mesh::SideSetEntry add_elem_side_pair(stk::mesh::Entity elem, int one_based_ordinal)
    {
        int zero_based_side_ordinal = convert_to_zero_based_ordinal(one_based_ordinal);
        return stk::mesh::SideSetEntry{elem, zero_based_side_ordinal};
    }

    void convert_elem_sides_pairs_into_sideset(const stk::mesh::BulkData& bulk, const std::vector<int>& elem_side, stk::mesh::SideSet& sideset)
    {
        for(size_t is=0; is<elem_side.size() / 2; ++is)
        {
            stk::mesh::Entity const elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elem_side[is*2]);
            if (bulk.is_valid(elem))
                sideset.push_back(add_elem_side_pair(elem, elem_side[is*2+1]));
        }
    }

    void convert_block_to_sideset(const stk::mesh::BulkData& bulk, const Ioss::SideBlock* block, stk::mesh::SideSet& sideset)
    {
        if (stk::io::include_entity(block))
        {
            std::vector<int> elem_side ;
            block->get_field_data("element_side", elem_side);
            convert_elem_sides_pairs_into_sideset(bulk, elem_side, sideset);
        }
    }

    void convert_ioss_sideset_to_stk_sideset(const stk::mesh::BulkData& bulk, const Ioss::SideSet* sset, stk::mesh::SideSet& sideset)
    {
        if(stk::io::include_entity(sset))
        {
            for (size_t i=0; i < sset->block_count(); i++)
            {
                Ioss::SideBlock *block = sset->get_block(i);
                convert_block_to_sideset(bulk, block, sideset);
            }
        }
    }

    void add_ioss_sideset_to_stk_sideset_with_sideset_id(const stk::mesh::BulkData& bulk, const Ioss::SideSet* sset, std::vector<IdAndSideSet>& sideset_data)
    {
        int sset_id = sset->get_property("id").get_int();
        sideset_data.push_back(IdAndSideSet{sset_id, stk::mesh::SideSet{}});
        stk::mesh::SideSet& currentSideset = sideset_data.back().sideSet;
        convert_ioss_sideset_to_stk_sideset(bulk, sset, currentSideset);
    }

public:
    void fill_sideset_data(std::vector<IdAndSideSet>& sidesets) const
    {
        sidesets = sideset_data;
    }

private:
    void extract_sideset_data_from_io(std::vector<IdAndSideSet> &sidesetData)
    {
        stk::mesh::BulkData &bulk = this->bulk_data();
        Ioss::Region *region = m_input_files[m_active_mesh_index]->get_input_io_region().get();
        for ( const Ioss::SideSet * sset : region->get_sidesets() )
            add_ioss_sideset_to_stk_sideset_with_sideset_id(bulk, sset, sidesetData);
    }

    std::vector<IdAndSideSet> sideset_data;
};

class BulkDataTester : public stk::mesh::BulkData
{
public:
    BulkDataTester(stk::mesh::MetaData & mesh_meta_data
                   , stk::ParallelMachine parallel)
    :   stk::mesh::BulkData(mesh_meta_data, parallel)
    {
    }

    void tester_save_sideset_data(int sideset_id, const stk::mesh::SideSet& data)
    {
        save_sideset_data(sideset_id, data);
    }
};

void write_exo_file(BulkDataTester &bulkData, const std::string &filename, const SideSetData& sidesetData);
void load_mesh_and_fill_sideset_data(StkMeshIoBrokerTester &stkIo, SideSetData &sidesetData);
void setup_io_broker_for_read(stk::io::StkMeshIoBroker &stkIo, stk::mesh::BulkData &bulkData, std::string filename, ReadMode readMode);
SideSetData get_sideset_data_from_written_file(BulkDataTester &bulk, const std::string& outputFileName);
void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& inputFileName, const std::string& outputFileName);
void compare_sidesets(const std::string& inputFileName,
                      BulkDataTester &bulk1,
                      const SideSetData &sidesetData1,
                      BulkDataTester &bulk2,
                      const SideSetData &sidesetData2);

}
}
}
#endif
