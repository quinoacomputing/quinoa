#include "mesh_adapter.h"

#include "Base/Exception.h"

namespace AMR {

    void mesh_adapter_t::init_node_store(coord_type* m_x, coord_type* m_y, coord_type* m_z, size_t* graph_size)
    {
        node_store.set_x(m_x);
        node_store.set_y(m_y);
        node_store.set_z(m_z);
        node_store.m_graphsize = graph_size;
    }

    void mesh_adapter_t::init_with_nodes(coord_type* m_x, coord_type* m_y, coord_type* m_z, size_t* graph_size)
    {
        init_node_store(m_x, m_y, m_z, graph_size);
        //init(); // TODO: This also needs to call init if you want node support to work
    }
    AMR::refinement_t
    mesh_adapter_t::init(const std::vector<size_t>& tetinpoel, size_t num_nodes) {
        node_connectivity.fill_initial_nodes(num_nodes);

        auto ref = AMR::refinement_t(tet_store, node_connectivity);

        consume_tets(tetinpoel);
        tet_store.generate_edges();
        return ref;
    }

    // This would be a candidate for a nice design pattern with runtime
    // selectable functionality..
    // TODO: Document this
    void mesh_adapter_t::consume_tets(const std::vector< std::size_t >& tetinpoel )
    {
        for (size_t i = 0; i < tetinpoel.size(); i+=4)
        {
            tet_t t = {
                {
                    tetinpoel[i],
                    tetinpoel[i+1],
                    tetinpoel[i+2],
                    tetinpoel[i+3]
                }
            };

            tet_store.add(t, AMR::Refinement_Case::initial_grid);
        }
    }

    /**
     * @brief Place holder function to evaluate error estimate at
     * nodes, and therefor mark things as needing to be refined
     */
    void mesh_adapter_t::evaluate_error_estimate() {
        // TODO: This could be implicit based of refinement_criteria value?
        for (auto& kv : tet_store.edge_store.edges)
        {
            // Mark them as needing refinement
            if (kv.second.refinement_criteria > refinement_cut_off)
            {
                kv.second.needs_refining = true;
            }
            else
            {
                // TODO: Check this won't be overwriting valuable
                // information from last iteration
                kv.second.needs_refining = false;
            }
        }
    }

    // TODO: Document this
    void mesh_adapter_t::uniform_refinement()
    {
        for (auto& kv : tet_store.edge_store.edges) {
            kv.second.refinement_criteria = 1.0;
        }
        evaluate_error_estimate();
        mark_refinement();
        perform_refinement();
    }

    /**
     * @brief Function to detect the compatibility class (1,
     * 2, or 3) based on the number of locked edges and the existence
     * of intermediate edges
     *
     * @param num_locked_edges The number of locked edges
     * @param refinement_case The refinement case of the tet
     *
     * @return The compatibility class of the current scenario
     */
    int mesh_adapter_t::detect_compatibility(int num_locked_edges,
            AMR::Refinement_Case refinement_case)
    {
        int compatibility = 0;

        // Only 1:2 and 1:4 are intermediates and eligible for class3
        if (
                (refinement_case == AMR::Refinement_Case::one_to_two) or
                (refinement_case == AMR::Refinement_Case::one_to_four)
           )
        {
            compatibility = 3;
        }
        else {
            if (num_locked_edges == 0) {
                compatibility = 1;
            }
            else {
                compatibility = 2;
            }
        }

        return compatibility;
    }

    /**
     * @brief Function which implements the main refinement algorithm from
     * the paper Iterating over the cells, deciding which refinement and
     * compatibility types are appropriate etc
     */
    void mesh_adapter_t::mark_refinement() {

        // Paper says the average actual num rounds will be 5-15
        const size_t max_num_rounds = 50;

        // Mark refinements
        size_t iter;
        //Iterate until convergence
        for (iter = 0; iter < max_num_rounds; iter++)
        {

            tet_store.marked_refinements.set_state_changed(false);

            // Loop over Tets.
            for (const auto& kv : tet_store.tets)
            {
                size_t tet_id = kv.first;

                //Only apply checks to tets on the active list
                if (tet_store.is_active(tet_id)) {
                    int compatibility = 1;
                    int num_locked_edges = 0;

                    // Loop over nodes and count the number which need refining
                    int num_to_refine = 0;

                    // This is useful for later inspection
                    edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);

                    // TODO: not all of this will be common to all three cases,
                    // so can be pushed down to the relevant case

                    //Iterate over edges
                    for(auto & key : edge_list)
                    {
                        //Count locked edges and edges in need of
                        // refinement
                        // Count Locked Edges
                        if(tet_store.edge_store.get(key).lock_case != AMR::Edge_Lock_Case::unlocked)
                        {
                            num_locked_edges++;
                        }
                        else
                        {
                            // Count edges which need refining
                            //  We check in here as we won't refine a
                            //  locked edge and will thus ignore it
                            if (tet_store.edge_store.get(key).needs_refining)
                            {
                                num_to_refine++;
                            }
                        }
                    }

                    AMR::Refinement_Case refinement_case = tet_store.get_refinement_case(tet_id);

                    //If we have any tets to refine
                    if (num_to_refine > 0)
                    {
                        //Determine compatibility case
                        compatibility = detect_compatibility(num_locked_edges,
                                refinement_case);

                        // Now check num_to_refine against situations
                        if (compatibility == 1)
                        {
                            refinement_class_one(num_to_refine, tet_id);
                        }
                        else if (compatibility == 2)
                        {
                            refinement_class_two(edge_list, tet_id);
                        }
                        else { // compatibility == 3
                            refinement_class_three(tet_id);
                        }

                        /*
                        // Write temp mesh out
                        std::string temp_file =  "temp." +
                        std::to_string(iter) + "." +
                        std::to_string(tet_id) + ".exo";

                        std::cout << "Writing " << temp_file << std::endl;
                        Adaptive_UnsMesh outmesh(
                        get_active_inpoel(), x(), y(), z()
                        );
                        tk::ExodusIIMeshWriter( temp_file, tk::ExoWriter::CREATE ).
                        writeMesh(outmesh);
                        */

                    } // if num_to_refine
                    else {
                        // If we got here, we don't want to refine this guy
                        tet_store.marked_refinements.add(tet_id, AMR::Refinement_Case::none);
                    }
                } // if active
            } // For

            // If nothing changed during that round, break
            if (!tet_store.marked_refinements.get_state_changed())
            {
                break;
            }
        }

    }

    // TODO: Document
    void mesh_adapter_t::print_tets() {
        tet_store.print_tets();
    }

    /**
     * @brief Function to call refinement after each tet has had it's
     * refinement case marked and calculated
     */
    void mesh_adapter_t::perform_refinement()
    {
        // Track tets which needs to be deleted this iteration
        std::set<size_t> round_two;

        // Do refinements
        for (const auto& kv : tet_store.tets)
        {
            size_t tet_id = kv.first;
            size_t parent_id = 0;

            if (tet_store.has_refinement_decision(tet_id))
            {
                switch(tet_store.marked_refinements.get(tet_id))
                {
                    case AMR::Refinement_Case::one_to_two:
                        refiner.refine_one_to_two(tet_id);
                        break;
                    case AMR::Refinement_Case::one_to_four:
                        refiner.refine_one_to_four(tet_id);
                        break;
                    case AMR::Refinement_Case::one_to_eight:
                        refiner.refine_one_to_eight(tet_id);
                        break;
                    case AMR::Refinement_Case::two_to_eight:
                        parent_id = tet_store.get_parent_id(tet_id);
                        round_two.insert(parent_id);
                        break;
                    case AMR::Refinement_Case::four_to_eight:
                        parent_id = tet_store.get_parent_id(tet_id);
                        round_two.insert(parent_id);
                        break;
                    case AMR::Refinement_Case::initial_grid:
                        // Do nothing
                    case AMR::Refinement_Case::none:
                        // Do nothing
                        break;
                        // No need for default as enum is explicitly covered?
                    //default:
                        //Assert(0);
                        //break;
                }
                // Mark tet as not needing refinement
                tet_store.marked_refinements.erase(tet_id);
            }
        }

        for (const auto i : round_two)
        {
            AMR::Refinement_State& element = tet_store.data(i);

            if (element.num_children == 2)
            {
                refiner.derefine_two_to_one(i);
            }
            else if (element.num_children == 4)
            {
                refiner.derefine_four_to_one(i);
            }
            else {
                std::cout << "num children " << element.num_children << std::endl;
                Assert(0, "Invalid number of children");
            }

            refiner.refine_one_to_eight(i);
            tet_store.unset_marked_children(i); // FIXME: This will not work well in parallel
            element.refinement_case = AMR::Refinement_Case::one_to_eight;
        }
        /*
           for (size_t t : tet_delete_list)
           {
        // TODO: also decrease number of children
        tets.erase(t);
        }
        */

        // Clean up dead edges
        // clean_up_dead_edges(); // Nothing get's marked as "dead" atm?

        std::cout << "Total Edges : " << tet_store.edge_store.size() << std::endl;
        std::cout << "Total Tets : " << tet_store.size() << std::endl;
        //std::cout << "Total Nodes : " << m_x.size() << std::endl;

        tet_store.print_node_types();
        //node_connectivity.print();
    }

    /**
     * @brief A method implementing "Algorithm 1" from the paper
     *
     * @param num_to_refine Number of edges to refine
     * @param tet_id The id of the given tet
     */
    void mesh_adapter_t::refinement_class_one(int num_to_refine, size_t tet_id)
    {
        // "If nrefine = 1
        // Accept as a 1:2 refinement"
        if (num_to_refine == 1)
        {
            //node_pair_t nodes = find_single_refinement_nodes(edge_list);
            //refine_one_to_two( tet_id, nodes[0], nodes[1]);
            // TODO: I'm pretty sure this needs to mark the edge for reifnement
            tet_store.mark_one_to_two(tet_id);
        }

        // "Else if nrefine = 2 OR nrefine = 3"
        else if (num_to_refine > 1 && num_to_refine < 4)
        {

            // We need to detect if the edges which need to refine are
            // on the same face
            // and if so which face so we know how to 1:4

            face_list_t face_list = tet_store.generate_face_lists(tet_id);
            bool edges_on_same_face = false;
            size_t face_refine_id = 0;

            // Iterate over each face
            for (size_t face = 0; face < NUM_TET_FACES; face++)
            {
                int num_face_refine_edges = 0;
                face_ids_t face_ids = face_list[face];

                edge_list_t face_edge_list = AMR::edge_store_t::generate_keys_from_face_ids(face_ids);
                // For this face list, see which ones need refining
                for (size_t k = 0; k < NUM_FACE_NODES; k++)
                {
                    edge_t key = face_edge_list[k];
                    if (tet_store.edge_store.get(key).needs_refining == true)
                    {
                        num_face_refine_edges++;
                    }
                }
                if (num_face_refine_edges == num_to_refine)
                {
                    edges_on_same_face = true;
                    face_refine_id = face;
                    break;
                }
            }

            // "If active edges are on the same face
            // Activate any inactive edges of the face
            // Accept as a 1:4 // refinement"
            if (edges_on_same_face)
            {
                //size_t opposite_offset = AMR::node_connectivity_t::face_list_opposite(face_list,
                        //face_refine_id);

                //tet_t tet = tet_store.get(tet_id);
                //size_t opposite_id = tet[opposite_offset];

                // Activate edges on this face
                edge_list_t face_edge_list = AMR::edge_store_t::generate_keys_from_face_ids(face_list[face_refine_id]);

                for (size_t k = 0; k < NUM_FACE_NODES; k++)
                {
                    edge_t key = face_edge_list[k];
                    tet_store.edge_store.mark_for_refinement(key);
                }

                //refiner.refine_one_to_four(tet_id, face_list[face_refine_id],
                //opposite_id);
                tet_store.mark_one_to_four(tet_id);
            }
            // "Else if active edges are not on the same face
            // Activate all edges
            // Accept as a 1:8 refinement"
            else {
                //refiner.refine_one_to_eight(tet_id);
                tet_store.mark_edges_for_refinement(tet_id);
                tet_store.mark_one_to_eight(tet_id);
            }

        }

        // "Else if nrefine > 3
        // Activate any inactive edges
        // Accept as a 1:8 refinement"
        else if (num_to_refine > 3)
        {
            //refiner.refine_one_to_eight(tet_id);
            tet_store.mark_one_to_eight(tet_id);
        }
    }

    // TODO: Document this
    void mesh_adapter_t::lock_tet_edges(size_t tet_id) {
        edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
        for (size_t k = 0; k < NUM_TET_EDGES; k++)
        {
            edge_t key = edge_list[k];
            if (tet_store.edge_store.get(key).lock_case == AMR::Edge_Lock_Case::unlocked)
            {
                tet_store.edge_store.get(key).lock_case = AMR::Edge_Lock_Case::locked;
            }
        }
    }

    // TODO: Document this
    // TODO: This has too similar a name to deactivate_tet
    void mesh_adapter_t::deactivate_tet_edges(size_t tet_id) {
        edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
        for (size_t k = 0; k < NUM_TET_EDGES; k++)
        {
            edge_t key = edge_list[k];

            tet_store.edge_store.unmark_for_refinement(key);
            // TODO: Should this also set something to not need derefining?
            tet_store.edge_store.get(key).needs_derefining = false;
        }
    }

    /**
     * @brief An implementation of "Algorithm 2" from the paper
     *
     * @param edge_list The list of edges for the given tet
     * @param tet_id The id of the given tet
     */
    void mesh_adapter_t::refinement_class_two(edge_list_t edge_list, size_t tet_id)
    {
        // "Deactivate all locked edges"

        // count number of active edges
        int num_active_edges = 0;
        for (size_t k = 0; k < NUM_TET_EDGES; k++)
        {
            edge_t key = edge_list[k];
            if (tet_store.edge_store.get(key).lock_case != AMR::Edge_Lock_Case::unlocked)
            {
                tet_store.edge_store.unmark_for_refinement(key);
            }
            // "Count number of active edges"
            if (tet_store.edge_store.get(key).needs_refining == true) {
                num_active_edges++;
            }
        }

        // Find out of two active edges live on the same face
        bool face_refine = false;
        size_t face_refine_id = 0; // FIXME: Does this need a better default
        face_list_t face_list = tet_store.generate_face_lists(tet_id);

        // Iterate over each face
        for (size_t face = 0; face < NUM_TET_FACES; face++)
        {
            int num_face_refine_edges = 0;

            face_ids_t face_ids = face_list[face];
            edge_list_t face_edge_list = AMR::edge_store_t::generate_keys_from_face_ids(face_ids);
            // For this face list, see which ones need refining
            for (size_t k = 0; k < NUM_FACE_NODES; k++)
            {
                edge_t key = face_edge_list[k];
                if (tet_store.edge_store.get(key).needs_refining == true)
                {
                    num_face_refine_edges++;
                }

                // Check for locked edges
                // This case only cares about faces with no locks
                if (tet_store.edge_store.get(key).lock_case != AMR::Edge_Lock_Case::unlocked)
                {
                    // Abort this face
                    num_face_refine_edges = 0;
                    continue;
                }
            }
            if (num_face_refine_edges >= 2)
            {
                face_refine = true;
                face_refine_id = face;
                break;
            }
        }

        // "If nrefine = 1
        // Accept as 1:2 refinement"
        if (num_active_edges == 1)
        {
            //node_pair_t nodes = find_single_refinement_nodes(edge_list);
            //refine_one_to_two( tet_id, nodes[0], nodes[1]);
            tet_store.mark_one_to_two(tet_id);
        }
        // "Else if any face has nrefine >= 2 AND no locked edges
        // Active any inactive edges of the face
        // Accept as a 1:4 refinement"
        else if (face_refine)
        {
            edge_list_t face_edge_list =
                AMR::edge_store_t::generate_keys_from_face_ids(face_list[face_refine_id]);

            for (size_t k = 0; k < NUM_FACE_NODES; k++)
            {
                edge_t key = face_edge_list[k];
                tet_store.edge_store.mark_for_refinement(key);
            }

            //refiner.refine_one_to_four(tet_id, face_list[face_refine_id], opposite_id);
            tet_store.mark_one_to_four(tet_id);
        }

        // "Else
        // Deactivate all edges
        // Mark all edges as locked"
        else {
            deactivate_tet_edges(tet_id);
            lock_tet_edges(tet_id);
        }

    }

    /**
     * @brief Based on a tet_id, decide if it's current state of locked
     * and marked edges maps to a valid refinement case. The logic for
     * this was derived from talking to JW and reading Chicoma.
     *
     * It basically just checks if something a 1:2 and has 3
     * intermediates and 3 makred edges, or is a 1:4 and has 5/6
     * intermediates
     *
     * @param child_id the id of the tet to check
     *
     * @return A bool saying if the tet is in a valid state to be refined
     */
    bool mesh_adapter_t::check_valid_refinement_case(size_t child_id) {

        edge_list_t edge_list = tet_store.generate_edge_keys(child_id);

        size_t num_to_refine = 0;
        size_t num_intermediate = 0;
        size_t unlocked = 0;
        size_t locked = 0;

        for (size_t k = 0; k < NUM_TET_EDGES; k++)
        {
            edge_t key = edge_list[k];

            // Count intermediate edges
            if (tet_store.edge_store.get(key).lock_case == AMR::Edge_Lock_Case::intermediate)
            {
                num_intermediate++;
            }

            // Count number of marked for refinement
            if (tet_store.edge_store.get(key).needs_refining)
            {
                num_to_refine++;
            }

            if (tet_store.edge_store.get(key).lock_case == AMR::Edge_Lock_Case::unlocked)
            {
                unlocked++;
            }

            if (tet_store.edge_store.get(key).lock_case == AMR::Edge_Lock_Case::locked)
            {
                locked++;
            }

        }

        AMR::Refinement_State& element = tet_store.data(child_id);

        // check if element is 1:2
        if (element.refinement_case == AMR::Refinement_Case::one_to_two)
        {
            // If so check it has 3 intermediates and 3 which need refining
            if (num_intermediate != 3 || num_to_refine != 3) {
                return false;
            }
        }

        // check if element is 1:4
        else if (element.refinement_case == AMR::Refinement_Case::one_to_four)
        {
            // TODO: Check if it's a center tet for a 1:4
            // FIXME: Is this even needed? How else would you get these
            // combinations? Can't we just combine these two checks?

            bool is_center_tet = tet_store.is_center(child_id);

            if (is_center_tet)
            {
                if (num_to_refine != 0 || num_intermediate != 6)
                {
                    return false;
                }
            }
            else { // Is one of the outsides (not center)
                if (num_to_refine != 1 || num_intermediate != 5)
                {
                    return false;
                }
            }
        }

        // If it makes it here, it's compatible
        return true;

    }

    /**
     * @brief Place holder method for the implementation of "Algorithm
     * 3" from the paper
     */
    // TODO: Does this parse a childs siblings multiple times?
    void mesh_adapter_t::refinement_class_three(size_t tet_id)
    {

        // "Identify parent element iparent"
        // TODO: WE should either always use the id to fetch, or always do the data lookup
        //size_t parent_id = master_elements.get_parent(tet_id);
        size_t parent_id = tet_store.get_parent_id(tet_id);

        // NOTE: This implies comms when we use these ids?
        child_id_list_t children = tet_store.data(parent_id).children;

        // "Do for each child element ielement
        // Activate all non-locked edges
        // Deactivate all locked edges"
        for (size_t i = 0; i < children.size(); i++)
        {
            // TODO: Is this in element or tet ids?
            edge_list_t edge_list = tet_store.generate_edge_keys(children[i]);
            for (size_t k = 0; k < NUM_TET_EDGES; k++)
            {
                edge_t key = edge_list[k];
                if (tet_store.edge_store.get(key).lock_case == AMR::Edge_Lock_Case::unlocked)
                {
                    tet_store.edge_store.mark_for_refinement(key);
                }
                else {
                    tet_store.edge_store.unmark_for_refinement(key);
                }
            }
        }

        // "Set compatible = TRUE
        bool compatible = true;

        // Do for each child element ielement
        // If ielement is not a valid refinement case
        // compatible = FALSE"
        for (size_t i = 0; i < children.size(); i++)
        {
            size_t child = children[i];
            if ( !check_valid_refinement_case(child) )
            {
                compatible = false;
            }
        }

        // "If compatible = FALSE
        // Do for each child element ielement
        // Deactive all edges of ielement
        // Mark all edges of ielement as locked
        // Mark ielement as normal"
        if (compatible == false)
        {
            for (size_t i = 0; i < children.size(); i++)
            {
                size_t child = children[i];
                deactivate_tet_edges(child);
                lock_tet_edges(child);
            }
        }
        else {
            // Accept as 2:8 or 4:8
            AMR::Refinement_State& element = tet_store.data(tet_id);
            if (element.refinement_case == AMR::Refinement_Case::one_to_two)
            {
                tet_store.mark_two_to_eight(tet_id);
            }
            else if (element.refinement_case == AMR::Refinement_Case::one_to_four)
            {
                tet_store.mark_four_to_eight(tet_id);
            }

        }
    }
}

