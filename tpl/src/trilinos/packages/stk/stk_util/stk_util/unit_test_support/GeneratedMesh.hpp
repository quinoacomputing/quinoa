/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_IO_UTIL_GENERATEDMESH_HPP
#define STK_IO_UTIL_GENERATEDMESH_HPP
#include <string>
#include <vector>

namespace stk {
  namespace io {
    namespace util {
      class GeneratedMesh
      {
      public:
	enum ShellLocation {MX =  0,  PX =  1,
			    MY =  2,  PY =  3,
			    MZ =  4,  PZ =  5};

	/**
	   Generate a cube mesh of size 'num_x' by 'num_y' by 'num_z' elements.
	   By default, the mesh is generated on a single processor.  If 'proc_count' is
	   greater than 1, then the mesh will be distributed over 'proc_count' processors
	   and this process will get the portion of the mesh for 'my_proc'.
	   The mesh will be decomposed along the 'Z' axis so 'num_z' must be greater than
	   or equal to 'proc_count' and for even distribution of the hexes 'num_z' mod 'proc_count'
	   should be zero.
	 
	   The mesh can optionally include shell elements along each face of the cube mesh.
	   These are specified via the 'add_shell_block' function.
	 
	   The mesh can optionally include nodesets/sidesets along each
	   face of the cube mesh.  These are specified via the
	   'add_nodesets' and 'add_sidesets' functions.
	 
	   If the 'parameters' string constructor is used, the string
	   is parsed to determine the intervals in each direction and,
	   optionally, additional information.  The form of the string
	   is "IxJxK" where I, J, and K are  the number of intervals
	   in the X, Y, and Z directions respectively and the "x" are
	   literal 'x' characters.  For example, the constructor
	   GeneratedMesh("10x12x14") will create the same mesh as
	   GeneratedMesh(10,12,14)

	   Additional valid options are:
           - help -- no argument, shows valid options
           - show -- no argument, prints out a summary of the
	   GeneratedMesh() parameters. The output will look similar
	   to:
	   \code
	   "10x12x8|shell:xX|bbox:-10,-10,-10,10,10,10|nodeset:xyz|sideset:XYZ|show"

	   Mesh Parameters:
	   Intervals: 10 by 12 by 8
	   X = 2       * (0..10) + -10     Range: -10 <= X <= 10
	   Y = 1.66667 * (0..12) + -10     Range: -10 <= Y <= 10
	   Z = 2.5     * (0..8)  + -10     Range: -10 <= Z <= 10
	   Node Count (total)    = 1287
	   Element Count (total) = 1152
	   Block Count           = 3
	   Nodeset Count         = 3
	   Sideset Count         = 3
	   \endcode

	   - shell -- argument = xXyYzZ which specifies whether there is a shell
	   block at that location. 'x' is minimum x face, 'X' is maximum x face,
	   similarly for y and z.  Note that the argument string is a single
	   multicharacter string.  You can add multiple shell blocks to a face,
	   for example, shell:xxx would add three layered shell blocks on the
	   minimum x face.  An error is output if a non xXyYzZ character is
	   found, but execution continues.

	   - nodeset -- argument = xXyYzZ which specifies whether there is
	   a nodeset at that location. 'x' is minimum x face, 'X' is
	   maximum x face, similarly for y and z.  Note that the argument
	   string is a single multicharacter string.  You can add multiple
	   nodesets to a face, for example, nodeset:xxx would add three
	   nodesets on the minimum x face.  An error is output if a non
	   xXyYzZ character is found, but execution continues.

	   - sideset -- argument = xXyYzZ which specifies whether there is
	   a sideset at that location. 'x' is minimum x face, 'X' is
	   maximum x face, similarly for y and z.  Note that the argument
	   string is a single multicharacter string.  You can add multiple
	   sidesets to a face, for example, sideset:xxx would add three
	   sidesets on the minimum x face.  An error is output if a non
	   xXyYzZ character is found, but execution continues.  If there
	   is a shell block specified on that face, then the sideset will
	   be on the shell elements; else the sideset will be on the hex
	   elements.

	   - zdecomp -- argument = n0, n1, n2, ..., n#proc-1 which are the number
	   of intervals in the z direction for each processor in a pallel run.
	   If this option is specified, then the total number of intervals in the
	   z direction is the sum of the n0, n1, ... An interval count must be
	   specified for each processor.  If this option is not specified, then
	   the number of intervals on each processor in the z direction is
	   numZ/numProc with the extras added to the lower numbered processors.

	   - scale -- argument = xs, ys, zs which are the scale factors in the x,
	   y, and z directions. All three must be specified if this option is
	   present.

	   - offset -- argument = xoff, yoff, zoff which are the offsets in the
	   x, y, and z directions.  All three must be specified if this option
	   is present.

	   - bbox -- argument = xmin, ymin, zmin, xmax, ymax, zmax
	   which specify the lower left and upper right corners of
	   the bounding box for the generated mesh.  This will
	   calculate the scale and offset which will fit the mesh in
	   the specified box.  All calculations are based on the currently
	   active interval settings. If scale or offset or zdecomp
	   specified later in the option list, you may not get the
	   desired bounding box.

	   - rotate -- argument = axis,angle,axis,angle,...
	   where axis is 'x', 'y', or 'z' and angle is the rotation angle in
	   degrees. Multiple rotations are cumulative. The composite rotation
	   matrix is applied at the time the coordinates are retrieved after
	   scaling and offset are applied.

	   The unrotated coordinate of a node at grid location i,j,k is:
	   \code
	   x = x_scale * i + x_off,
	   y = z_scale * j + y_off,
	   z = z_scale * k + z_off,
	   \endcode

	   The extent of the unrotated mesh will be:
	   \code
	   x_off <= x <= x_scale * numX + x_off
	   y_off <= y <= y_scale * numY + y_off
	   z_off <= z <= z_scale * numZ + z_off
	   \endcode

	   If an unrecognized option is specified, an error message will be
	   output and execution will continue.

	   An example of valid input is:
	   \code
	   "10x20x40|scale:1,0.5,0.25|offset:-5,-5,-5|shell:xX"
	   \endcode

	   This would create a mesh with 10 intervals in x, 20 in y, 40 in z
	   The mesh would be centered on 0,0,0 with a range of 10 in each
	   direction. There would be a shell layer on the min and max
	   x faces.

	   NOTE: All options are processed in the order they appear in
	   the parameters string (except rotate which is applied at the
	   time the coordinates are generated/retrieved)
	*/
	explicit GeneratedMesh(const std::string &parameters, int proc_count = 1, int my_proc = 0);
	GeneratedMesh(int num_x, int num_y, int num_z,        int proc_count = 1, int my_proc = 0);
	~GeneratedMesh();

	/**
	 * Add a shell block along the specified face of the hex mesh.
	 * The shell blocks will maintain the order of definition. The
	 * first shell block defined will be block 2; the hex block has id
	 * 1.  The loc options are:
	 * - MX = add shell layer on the face with minimum X 
	 * - PX = add shell layer on the face with maximum X 
	 * - MY = add shell layer on the face with minimum Y 
	 * - PY = add shell layer on the face with maximum Y 
	 * - MZ = add shell layer on the face with minimum Z 
	 * - PZ = add shell layer on the face with maximum Z
	 *
	 */
	size_t add_shell_block(ShellLocation loc);
    
	/**
	 * Add a nodeset along the specified face of the hex mesh.
	 * The nodesets will maintain the order of definition. The
	 * first nodeset defined will be nodeset 1.
	 * The loc options are:
	 * - MX = add nodeset on the face with minimum X 
	 * - PX = add nodeset on the face with maximum X 
	 * - MY = add nodeset on the face with minimum Y 
	 * - PY = add nodeset on the face with maximum Y 
	 * - MZ = add nodeset on the face with minimum Z 
	 * - PZ = add nodeset on the face with maximum Z
	 *
	 */
	size_t add_nodeset(ShellLocation loc);
    
	/**
	 * Add a sideset along the specified face of the hex mesh.
	 * The sidesets will maintain the order of definition. The
	 * first sideset defined will be sideset 1. If there is a shell
	 * block specified on that face, then the sideset will be on the
	 * shell elements; otherwise the sideset will be on the hex
	 * elements. 
	 * The loc options are:
	 * - MX = add sideset on the face with minimum X 
	 * - PX = add sideset on the face with maximum X 
	 * - MY = add sideset on the face with minimum Y 
	 * - PY = add sideset on the face with maximum Y 
	 * - MZ = add sideset on the face with minimum Z 
	 * - PZ = add sideset on the face with maximum Z
	 *
	 */
	size_t add_sideset(ShellLocation loc);
    
	/**
	 * Specify the coordinate scaling and offset in all three
	 * spatial dimensions.
	 *
	 * node location of node at (i,j,k) is
	 * \code
	 * X = scale X * i + offset X
	 * Y = scale Y * i + offset Y
	 * Z = scale Z * i + offset Z
	 * \endcode
	 *
	 * WARNING: Should be called before retrieving node
	 * coordinates.
	 */
	void set_scale(double scl_x, double scl_y, double scl_z);
	void set_offset(double off_x, double off_y, double off_z);
	void set_bbox(double xmin, double ymin, double zmin,
		      double xmax, double ymax, double zmax);
	
	/**
	 * Set rotation.  Multiple calls are cumulative.
	 * Rotate 'angle_degrees' degrees about the axis 'axis'
	 * Center of rotation is about the origin and operates
	 * on the scaled/offset coordinates of the mesh.
	 */
	void set_rotation(const std::string &axis, double angle_degrees);

	/**
	 * Return number of nodes in the entire model.
	 */
	size_t node_count() const;

	/**
	 * Return number of nodes on this processor.
	 */
	size_t node_count_proc() const;

	/**
	 * Return number of element blocks in the entire model.
	 */
	size_t block_count() const;

	/**
	 * Return number of nodesets in the entire model.
	 */
	size_t nodeset_count() const;

	/**
	 * Return number of nodeset nodes on nodeset 'id'
	 */
	size_t nodeset_node_count(size_t id) const;

	/**
	 * Return number of nodeset nodes on nodeset 'id' on the current processor
	 */
	size_t nodeset_node_count_proc(size_t id) const;

	/**
	 * Return number of sidesets in the entire model.
	 */
	size_t sideset_count() const;

	/**
	 * Return number of sideset 'sides' on sideset 'id'
	 */
	size_t sideset_side_count(size_t id) const;

	/**
	 * Return number of sideset 'sides' on sideset 'id' on the current
	 * processor.
	 */
	size_t sideset_side_count_proc(size_t id) const;

	/**
	 * Return number of elements in all element blocks in the model.
	 */
	size_t element_count() const;

	/**
	 * Return number of shell elements in all element blocks in the model.
	 */
	size_t shell_element_count(ShellLocation) const;

	/**
	 * Return number of elements in all element blocks on this processor.
	 */
	size_t element_count_proc() const;

	/**
	 * Return number of shell elements in all element blocks on this processor.
	 */
	size_t shell_element_count_proc(ShellLocation) const;

	/**
	 * Return number of elements in the element block with id
	 * 'block_number'. The 'block_number' ranges from '1' to
	 * 'block_count()'. 
	 */
	size_t element_count(size_t block_number) const;

	/**
	 * Return number of elements on this processor in the element
	 * block with id 'block_number'. The 'block_number' ranges from
	 * '1' to 'block_count()'.
	 */
	size_t element_count_proc(size_t block_number) const;

	/**
	 * Returns pair containing "topology type string" and "number of
	 * nodes / element". The topology type string will be "hex8" for
	 * the hex element block and "shell4" for the shell element blocks.
	 */
	std::pair<std::string, int>  topology_type(size_t block_number) const;
    
	size_t communication_node_count_proc() const;
	void node_communication_map(std::vector<int> &map, std::vector<int> &proc);
    
	/** 
	 * Fill the passed in 'map' argument with the node map
	 * "map[local_position] = global_id" for the nodes on this
	 * processor.
	 */
	void node_map(std::vector<int> &map);

	/** 
	 * Fill the passed in 'map' argument with the element map
	 * "map[local_position] = global_id" for the elements on this
	 * processor in block "block_number".
	 */
	void element_map(size_t block_number, std::vector<int> &map) const;

	/** 
	 * Fill the passed in 'map' argument with the element map
	 * "map[local_position] = global_id" for all elements on this
	 * processor 
	 */
	void element_map(std::vector<int> &map) const;
    
	/** 
	 * Fill the passed in 'map' argument with the element map pair
	 * "map[local_position] = element global_id" and
	 * "map[local_position+1] = element local face id (0-based)" for
	 * all elements on the current processor having a face on the
	 * surface defined by ShellLocation.
	 */
	void element_surface_map(ShellLocation loc, std::vector<int> &map) const;
    
	/**
	 * Return the connectivity for the elements on this processor in
	 * the block with id 'block_number'. If the elements in this block
	 * have 'npe' nodes per element, then the first 'npe' entries in
	 * the 'conn' vector will be the nodal connectivity for the first
	 * element; the next 'npe' entries are the nodal connectivity for
	 * the second element.  The 'connect' vector will be resized to the
	 * size required to contain the nodal connectivity for the
	 * specified block; all information in 'connect' will be overwritten.
	 */
	void connectivity(size_t block_number, std::vector<int> &connect) const;

	/**
	 * Return the coordinates for all nodes on this processor.  The
	 * first 3 entries in the 'coord' vector are the x, y, and z
	 * coordinates of the first node, etc.  The 'coord' vector will be
	 * resized to the size required to contain the nodal coordinates;
	 * all information in 'coord' will be overwritten.
	 */
	void coordinates(std::vector<double> &coord) const;

	/**
	 * Return the coordinates for all nodes on this processor in
	 * separate vectors. The vectors will be resized to the size
	 * required to contain the nodal coordinates; all information in
	 * the vectors will be overwritten.
	 */
	void coordinates(std::vector<double> &x,
			 std::vector<double> &y,
			 std::vector<double> &z) const;
    
	/**
	 * Return the list of nodes in nodeset 'id' on this processor.
	 * The 'nodes' vector will be resized to the size required to
	 * contain the node list. The ids are global ids.
	 */
	void nodeset_nodes(size_t nset_id, std::vector<int> &nodes) const;

	/**
	 * Return the list of the face/ordinal pairs
	 * "elem_sides[local_position]   = element global_id" and
	 * "elem_sides[local_position+1] = element local face id (0-based)"
	 * for the faces in sideset 'id' on this
	 * processor.  The 'elem_sides' vector will be resized to the size
	 * required to contain the list. The element ids are global ids,
	 * the side ordinal is 0-based.
	 */
	void sideset_elem_sides(size_t nset_id, std::vector<int> &elem_sides) const;

      private:
    
	GeneratedMesh( const GeneratedMesh & );
	GeneratedMesh & operator = ( const GeneratedMesh & );

	void parse_options(const std::vector<std::string> &groups);
	void show_parameters() const;
	void initialize();
	
        std::vector<ShellLocation> shellBlocks;
	std::vector<ShellLocation> nodesets;
	std::vector<ShellLocation> sidesets;
        double rotmat[3][3];
	size_t numX, numY, numZ;
	size_t myNumZ, myStartZ;

	size_t processorCount;
	size_t myProcessor;

	double offX, offY, offZ; /** Offsets in X, Y, and Z directions */
	double sclX, sclY, sclZ; /** Scale in X, Y, and Z directions
				  * location of node at (i,j,k)
				  * position is (sclX*i+offX,
				  * sclY*i+offY, sclZ*i+offZ) */
	bool doRotation;
      };
    }
  }
}
#endif
