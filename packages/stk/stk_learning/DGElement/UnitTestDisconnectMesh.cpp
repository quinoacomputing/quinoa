#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>

#include "TetFixture.hpp"

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "../../stk_io/stk_io/StkMeshIoBroker.hpp"

class DisconnectMesh : public DGTetFixture {};

void disconnectMesh(const stk::mesh::MetaData &oldMeta,
                         stk::mesh::BulkData& oldBulkData,
                         stk::mesh::Selector subMeshSelector,
                         stk::mesh::MetaData &newMeta,
                         stk::mesh::BulkData &newBulkData);

TEST_F(DisconnectMesh, tet)
{
    std::vector<stk::mesh::EntityIdVector> tet_conn = {
            {1, 2, 3, 4}, // id 1
            {2, 3, 4, 5}  // id 2
    };

    std::vector< std::vector<double> > node_coords= {
            {0, 0, 0}, // 1
            {1, 0, 0}, // 2
            {0, 1, 0}, // 3
            {0.5, 0.5, 1.0}, // 6...just kidding, it's 4
            {1.0, 1.0, 1.0}
    };

    setup_mesh(tet_conn, node_coords);

    //////////////////////////////////////////////////////////////////////////////////////

    stk::mesh::MetaData newMetaData(get_meta().spatial_dimension());
    stk::mesh::BulkData newBulkData(newMetaData, get_bulk().parallel());
    disconnectMesh(get_meta(), get_bulk(), get_meta().universal_part(), newMetaData, newBulkData);

    stk::unit_test_util::write_mesh_using_stk_io("mike_new.g", newBulkData);
}

class TOSDTWD : public stk::unit_test_util::MeshFixture
{
protected:
    TOSDTWD() : stk::unit_test_util::MeshFixture(3) {}
    virtual ~TOSDTWD() {}
};

TEST_F(TOSDTWD, disconnect_mesh)
{
    std::string exodusFileName = unitTestUtils::getOption("-i", "generated:10x10x10");
    std::string mesh_name = unitTestUtils::getOption("-o", "blown_up.g");

    {
        setup_mesh(exodusFileName, stk::mesh::BulkData::NO_AUTO_AURA);
        stk::mesh::MetaData newMetaData(get_meta().spatial_dimension());
        stk::mesh::BulkData newBulkData(newMetaData, get_bulk().parallel());
        disconnectMesh(get_meta(), get_bulk(), get_meta().universal_part(), newMetaData, newBulkData);
        stk::unit_test_util::write_mesh_using_stk_io(mesh_name, newBulkData);
    }
}

TEST_F(TOSDTWD, birds)
{
    std::string velocityx = unitTestUtils::getOption("-vx", "10");
    std::string velocityy = unitTestUtils::getOption("-vy", "10");
    std::string steps = unitTestUtils::getOption("-steps", "100");
    std::string ax = unitTestUtils::getOption("-ax", "0");
    std::string ay = unitTestUtils::getOption("-ay", "9.8");
    std::string exodusFileName = unitTestUtils::getOption("-i", "generated:10x10x10");
    std::string results_name = "output_results.g";

    double velx = 10;
    {
        std::istringstream is(velocityx);
        is >> velx;
    }

    double vely = 10;
    {
        std::istringstream is(velocityy);
        is >> vely;
    }

    unsigned int num_steps;
    {
        std::istringstream is(steps);
        is >> num_steps;
    }

    double accelx = 0;
    {
        std::istringstream is(ax);
        is >> accelx;
    }

    double accely = 9.8;
    {
        std::istringstream is(ay);
        is >> accely;
    }

    std::cerr << "Using (vx,vy) and (gx,gy) of (" << velx << ", " << vely << ") and (" << ax << ", " << ay << ")" << std::endl;

    {
           stk::io::StkMeshIoBroker stkIo(get_comm());
           stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);

           stkIo.create_input_mesh();

           typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> CoordFieldType;

           CoordFieldType& field = stkIo.meta_data().declare_field<CoordFieldType>(stk::topology::NODE_RANK, "disp");
           CoordFieldType& velocity = stkIo.meta_data().declare_field<CoordFieldType>(stk::topology::ELEM_RANK, "velocity");
           stk::mesh::Field<double>& active_status = stkIo.meta_data().declare_field<stk::mesh::Field<double>>(stk::topology::ELEMENT_RANK, "active_status");

           double init_status = 1;
           std::vector<double> init_vec = {0,0,0};

           stk::mesh::put_field(field, stkIo.meta_data().universal_part(), init_vec.data());
           stk::mesh::put_field(velocity, stkIo.meta_data().universal_part(), init_vec.data());
           stk::mesh::put_field(active_status, stkIo.meta_data().universal_part(), &init_status);

           stkIo.populate_bulk_data();

           size_t fh = stkIo.create_output_mesh(results_name, stk::io::WRITE_RESULTS);

           stkIo.add_field(fh, field); /*@\label{io:results:add_field}*/
           stkIo.add_field(fh, velocity); /*@\label{io:results:add_field}*/
           stkIo.add_field(fh, active_status); /*@\label{io:results:add_field}*/

           std::vector<stk::mesh::Entity> elements;
           stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEM_RANK, elements);

           CoordFieldType *coords = stkIo.meta_data().get_field<CoordFieldType>(stk::topology::NODE_RANK, "coordinates");
           std::vector<std::vector<double> > direct(elements.size());

           double minx = 0, miny = 0, minz = 0, maxx = 0, maxy = 0, maxz = 0;

           double vx = velx, vy = vely, vz = 0;
           double gravity = accely;
           double gx = -0.5*accelx, gy = -gravity*0.5, gz = 0;
           double xc = 0, yc = 0, zc = 0;

           for(size_t i = 0; i < elements.size(); i++)
           {
               double *elem_vel = stk::mesh::field_data(velocity, elements[i]);

               elem_vel[0] = vx;
               elem_vel[1] = vy;
               elem_vel[2] = vz;

               unsigned num_nodes = stkIo.bulk_data().num_nodes(elements[i]);
               const stk::mesh::Entity *nodes = stkIo.bulk_data().begin_nodes(elements[i]);
               double exc = 0, eyc = 0, ezc = 0;
               for(unsigned j=0;j<num_nodes;++j)
               {
                   double *node_data = stk::mesh::field_data(*coords, nodes[j]);
                   exc += node_data[0];
                   eyc += node_data[1];
                   ezc += node_data[2];
                   minx = std::min(minx, node_data[0]);
                   miny = std::min(miny, node_data[1]);
                   minz = std::min(minz, node_data[2]);
                   maxx = std::max(maxx, node_data[0]);
                   maxy = std::max(maxy, node_data[1]);
                   maxz = std::max(maxz, node_data[2]);
               }

               exc /= num_nodes;
               eyc /= num_nodes;
               ezc /= num_nodes;
               xc += exc;
               yc += eyc;
               zc += ezc;
               direct[i].resize(3);
               direct[i][0] = exc;
               direct[i][1] = eyc;
               direct[i][2] = ezc;
           }

           xc /= elements.size();
           yc /= elements.size();
           zc /= elements.size();

           for(size_t i=0;i<direct.size();++i)
           {
               direct[i][0] -= xc;
               direct[i][1] -= yc;
               direct[i][2] -= zc;
               double mag = sqrt(direct[i][0]*direct[i][0] +
                                 direct[i][1]*direct[i][1] +
                                 direct[i][2]*direct[i][2]);
               direct[i][0] /= mag;
               direct[i][1] /= mag;
               direct[i][2] /= mag;
           }

           std::cerr << "(x,y,z) centroid = (" << xc << ", " << yc << ", " << zc << ")\n";

           double time_step = 0.1;

           for(unsigned int step = 0; step < num_steps; step++)
            {
                double time = step * time_step;

                for(size_t i = 0; i < elements.size(); i++)
                {
                    unsigned num_nodes = stkIo.bulk_data().num_nodes(elements[i]);
                    const stk::mesh::Entity *nodes = stkIo.bulk_data().begin_nodes(elements[i]);

                    bool are_any_nodes_in_contact = false;
                    bool left_half = false;
                    bool top_half = false;

                    std::vector<double> direction = {0,0,0};

                    double *status = stk::mesh::field_data(active_status, elements[i]);

                    double *elem_vel = stk::mesh::field_data(velocity, elements[i]);

                    elem_vel[0] += gx * time_step;
                    elem_vel[1] += gy * time_step;
                    elem_vel[2] += gz * time_step;

                    for(unsigned j=0;j<num_nodes;++j)
                    {
                        double *disp = stk::mesh::field_data(field, nodes[j]);

                        disp[0] += elem_vel[0] * time_step + gx * time_step*time_step;
                        disp[1] += elem_vel[1] * time_step + gy * time_step*time_step;
                        disp[2] += elem_vel[2] * time_step + gz * time_step*time_step;

                        double *node_data = stk::mesh::field_data(*coords, nodes[j]);
                        if(elem_vel[1]<0 && (node_data[1]+disp[1])<1e-6 && *status == 1)
                        {
                            are_any_nodes_in_contact = true;
                            if(node_data[0]<xc)
                                left_half = true;

                            if(node_data[1]>yc)
                                top_half = true;

                            direction[0] = node_data[0] - xc;
                            direction[1] = node_data[1] - xc;
                            direction[2] = node_data[2] - xc;
                            double mag = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
                            direction[0] /= mag;
                            direction[1] /= mag;
                            direction[2] /= mag;
                        }
                    }

                    if (are_any_nodes_in_contact)
                    {
                        double vel_mag = sqrt(elem_vel[0]*elem_vel[0] + elem_vel[1]*elem_vel[1]);

                        if(direction[1]<0)
                            direction[1] *=-1;

                        double coeff_restitution = 1.2;
                        if(top_half)
                            coeff_restitution = 1.1;

                        *status = 0;
                        elem_vel[0] = coeff_restitution*vel_mag*direction[0];
                        elem_vel[1] = coeff_restitution*vel_mag*direction[1];

                        if(left_half)
                            continue;
//                        elem_vel[1] = -coeff_restitution*elem_vel[1];
//                        if(left_half)
//                            elem_vel[0] = -coeff_restitution*elem_vel[0];
                    }
                    else
                    {
                        *status = 1;
                    }
                }

                stkIo.begin_output_step(fh, time);
                stkIo.write_defined_output_fields(fh);
                stkIo.end_output_step(fh);
            }
         }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void copyPartsToNewMesh(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::PartVector &allparts = oldMeta.get_mesh_parts();
    for(size_t i = 0; i < allparts.size(); i++)
    {
        stk::mesh::Part *part = NULL;
        if(allparts[i]->topology() != stk::topology::INVALID_TOPOLOGY)
        {
            part = &newMeta.declare_part_with_topology(allparts[i]->name(), allparts[i]->topology());
        }
        else
        {
            part = &newMeta.declare_part(allparts[i]->name());
        }
        if(stk::io::is_part_io_part(*allparts[i]))
        {
            stk::io::put_io_part_attribute(*part);
        }
    }
}

void copyFieldsToNewMesh(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::FieldVector &fields = oldMeta.get_fields();
    for(size_t i = 0; i < fields.size(); i++)
    {
        stk::mesh::FieldBase* newField = newMeta.declare_field_base(fields[i]->name(),
                                                                    fields[i]->entity_rank(),
                                                                    fields[i]->data_traits(),
                                                                    fields[i]->field_array_rank(),
                                                                    fields[i]->dimension_tags(),
                                                                    fields[i]->number_of_states());

        stk::mesh::Selector selectFieldParts = stk::mesh::selectField(*fields[i]);
        stk::mesh::PartVector oldParts;
        selectFieldParts.get_parts(oldParts);
        if(!oldParts.empty())
        {
            stk::mesh::PartVector newParts(oldParts.size());
            for(size_t k = 0; k < oldParts.size(); k++)
            {
                newParts[k] = &newMeta.get_part(oldParts[k]->mesh_meta_data_ordinal());
            }
            stk::mesh::Selector selectNewParts = stk::mesh::selectUnion(newParts);
            stk::mesh::put_field(*newField, selectNewParts, fields[i]->max_size(fields[i]->entity_rank()));
        }
    }
}

void copySelectedEntitiesToNewMesh(const stk::mesh::BulkData& oldBulkData,
                                   const stk::mesh::Selector subMeshSelector,
                                   std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                                   stk::mesh::MetaData &newMeta,
                                   stk::mesh::BulkData &newBulkData)
{
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = oldBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t j = 0; j < bucket.size(); j++)
            {
                const stk::mesh::PartVector &oldParts = bucket.supersets();
                stk::mesh::PartVector newParts(oldParts.size(), 0);
                for(size_t k = 0; k < oldParts.size(); k++)
                {
                    newParts[k] = &newMeta.get_part(oldParts[k]->mesh_meta_data_ordinal());
                }
                stk::mesh::Entity newEntity = newBulkData.declare_entity(rank, oldBulkData.identifier(bucket[j]), newParts);
                oldToNewEntityMap[bucket[j]] = newEntity;
            }
        }
    }
}

void copyRelationsToNewMesh(const stk::mesh::BulkData& oldBulkData,
                            const stk::mesh::Selector subMeshSelector,
                            const std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                            stk::mesh::BulkData &newBulkData)
{
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = oldBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t j = 0; j < bucket.size(); j++)
            {
                stk::mesh::Entity newFromEntity = oldToNewEntityMap.find(bucket[j])->second;
                for(stk::mesh::EntityRank downRank = stk::topology::NODE_RANK; downRank < rank; downRank++)
                {
                    unsigned numConnectedEntities = oldBulkData.num_connectivity(bucket[j], downRank);
                    const stk::mesh::Entity *connectedEntities = oldBulkData.begin(bucket[j], downRank);
                    for(unsigned connectOrder = 0; connectOrder < numConnectedEntities; connectOrder++)
                    {
                        stk::mesh::Entity newToEntity = oldToNewEntityMap.find(connectedEntities[connectOrder])->second;
                        newBulkData.declare_relation(newFromEntity, newToEntity, connectOrder);
                    }
                }
            }
        }
    }
}

void copyFieldDataToNewMesh(const stk::mesh::MetaData &oldMeta,
                            const stk::mesh::BulkData& oldBulkData,
                            const stk::mesh::Selector subMeshSelector,
                            const std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                            stk::mesh::MetaData &newMeta,
                            const stk::mesh::BulkData& newBulkData,
                            const stk::mesh::EntityVector &node_in_old_mesh)
{
    const stk::mesh::FieldVector &oldFields = oldMeta.get_fields();
    const stk::mesh::FieldVector &newFields = newMeta.get_fields();
    for(stk::mesh::EntityRank rank = stk::topology::EDGE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = oldBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t k = 0; k < oldFields.size(); k++)
            {
                if(bucket.field_data_is_allocated(*oldFields[k]))
                {
                    for(size_t j = 0; j < bucket.size(); j++)
                    {
                        stk::mesh::Entity oldEntity = bucket[j];
                        stk::mesh::Entity newEntity = oldToNewEntityMap.find(oldEntity)->second;
                        void *oldData = stk::mesh::field_data(*oldFields[k], oldEntity);
                        void *newData = stk::mesh::field_data(*newFields[k], newEntity);
                        memcpy(newData, oldData, stk::mesh::field_bytes_per_entity(*oldFields[k], oldEntity));
                    }
                }
            }
        }
    }

    const stk::mesh::BucketVector &buckets = newBulkData.get_buckets(stk::topology::NODE_RANK, subMeshSelector);
    for(size_t i = 0; i < buckets.size(); i++)
    {
        stk::mesh::Bucket &bucket = *buckets[i];
        for(size_t k = 0; k < newFields.size(); k++)
        {
            if(bucket.field_data_is_allocated(*newFields[k]))
            {
                for(size_t j = 0; j < bucket.size(); j++)
                {
                    stk::mesh::Entity newEntity = bucket[j];
                    stk::mesh::Entity oldEntity = node_in_old_mesh[newBulkData.identifier(newEntity)-1];
                    void *oldData = stk::mesh::field_data(*oldFields[k], oldEntity);
                    void *newData = stk::mesh::field_data(*newFields[k], newEntity);
                    memcpy(newData, oldData, stk::mesh::field_bytes_per_entity(*oldFields[k], oldEntity));
                }
            }
        }
    }
}


void disconnectMesh(const stk::mesh::MetaData &oldMeta,
                         stk::mesh::BulkData& oldBulkData,
                         stk::mesh::Selector subMeshSelector,
                         stk::mesh::MetaData &newMeta,
                         stk::mesh::BulkData &newBulkData)
{
    copyPartsToNewMesh(oldMeta, newMeta);
    copyFieldsToNewMesh(oldMeta, newMeta);

    newMeta.commit();

    const stk::mesh::PartVector &allOldparts = oldMeta.get_mesh_parts();
    const stk::mesh::PartVector &allNewparts = newMeta.get_mesh_parts();
    ThrowRequireMsg(allOldparts.size() == allNewparts.size(), "Error. Contact support.");

    std::map<stk::mesh::Entity, stk::mesh::Entity> oldToNewEntityMap;

    stk::mesh::EntityVector node_in_old_mesh;
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(oldMeta.locally_owned_part(), oldBulkData.buckets(stk::topology::ELEM_RANK), elements);
    size_t num_nodes_needed = 0;
    for(stk::mesh::Entity element : elements )
    {
        unsigned num_nodes_this_element = oldBulkData.num_nodes(element);
        const stk::mesh::Entity* nodes = oldBulkData.begin_nodes(element);
        num_nodes_needed += num_nodes_this_element;
        for(unsigned int i=0;i<num_nodes_this_element;++i)
        {
            node_in_old_mesh.push_back(nodes[i]);
        }
    }

    std::vector<size_t> num_nodes_needed_per_proc(oldBulkData.parallel_size(),0);
    num_nodes_needed_per_proc[oldBulkData.parallel_rank()] = num_nodes_needed;
    stk::all_reduce_sum(oldBulkData.parallel(), num_nodes_needed_per_proc.data(), num_nodes_needed_per_proc.data(), num_nodes_needed_per_proc.size());

    size_t offset_this_proc = 0;
    for(int i=0;i<oldBulkData.parallel_rank();++i)
        offset_this_proc += num_nodes_needed_per_proc[i];

    stk::mesh::EntityVector new_nodes(num_nodes_needed);
    newBulkData.modification_begin();
    for(size_t i=0;i<num_nodes_needed;++i)
    {
        const stk::mesh::PartVector parts = oldBulkData.bucket(node_in_old_mesh[i]).supersets();
        new_nodes[i] = newBulkData.declare_entity(stk::topology::NODE_RANK, i+offset_this_proc+1, parts);
    }
    newBulkData.modification_end();

    newBulkData.modification_begin();
    size_t node_counter = 0;
    for(stk::mesh::Entity element : elements)
    {
        stk::mesh::PartVector parts = oldBulkData.bucket(element).supersets();

        int index_root_part = -1;
        for(size_t i=0;i<parts.size();++i)
        {
            if(stk::mesh::is_topology_root_part(*parts[i]))
            {
                index_root_part = i;
                break;
            }
        }
        ThrowRequireMsg(index_root_part !=-1, "Oops");

        std::swap(parts[0], parts[index_root_part]);

        stk::mesh::EntityIdVector node_ids(oldBulkData.num_nodes(element));
        for(size_t i=0;i<node_ids.size();++i)
            node_ids[i] = node_counter + i + offset_this_proc + 1;
        node_counter += node_ids.size();
        stk::mesh::declare_element(newBulkData, parts, oldBulkData.identifier(element), node_ids);
    }
    newBulkData.modification_end();

    copyFieldDataToNewMesh(oldMeta, oldBulkData, subMeshSelector, oldToNewEntityMap, newMeta, newBulkData, node_in_old_mesh);
}
