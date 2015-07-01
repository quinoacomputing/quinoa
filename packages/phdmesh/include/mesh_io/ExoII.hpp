/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#ifndef phdmesh_ExoII_hpp
#define phdmesh_ExoII_hpp

#include <mesh/Types.hpp>
#include <mesh/FieldTraits.hpp>

namespace phdmesh {
namespace exodus {

class FilePart ;

//----------------------------------------------------------------------

struct ElementAttributes : public ArrayDimTag {
  const char * name() const ;
  static const ElementAttributes & tag();
private:
  ElementAttributes() {}
  ElementAttributes( const ElementAttributes & );
  ElementAttributes & operator = ( const ElementAttributes & );
};

struct GlobalLocalIndex : public ArrayDimTag {
  const char * name() const ;

  std::string to_string( unsigned size , unsigned index ) const ;

  unsigned to_index( unsigned size , const std::string & ) const ;

  static const GlobalLocalIndex & tag();
private:
  GlobalLocalIndex() {}
  GlobalLocalIndex( const GlobalLocalIndex & );
  GlobalLocalIndex & operator = ( const GlobalLocalIndex & );
};

class FileSchema {
public:

  typedef Field<double,ElementAttributes> AttributeField ;
  typedef Field<int,GlobalLocalIndex>     IndexField ;

  ~FileSchema();

  FileSchema( MetaData          & arg_schema ,
              const FieldBase       & arg_node_coordinates ,
              const AttributeField  & arg_elem_attributes ,
              const unsigned          arg_writer_rank = 0 );

  FileSchema( MetaData          & arg_schema ,
              const FieldBase       & arg_node_coordinates ,
              const AttributeField  & arg_elem_attributes ,
              const std::string     & arg_file_path ,
              ParallelMachine         arg_comm ,
              const unsigned          arg_writer_rank = 0 );

  /** Declare element part, default number of attributes. */
  void declare_part( Part & arg_part , int arg_id );

  /** Declare a node, edge, or face part */
  void declare_part( Part     & arg_part ,
                     int        arg_id ,
                     EntityType arg_type );

  /** Assign contiguous global indices [1..#] to nodes and elements.
   *  Elements are ordered by element block and then by identifier.
   */
  void assign_indices( BulkData & ) const ;

  MetaData          & m_schema ;
  const unsigned          m_io_rank ;
  const unsigned          m_dimension ;
  const FieldBase       & m_field_node_coord ;
  const AttributeField  & m_field_elem_attr ;
  const IndexField      & m_field_index ;

  const std::vector<const FilePart*> & parts( EntityType t ) const
    { return m_parts[t] ; }

private:
  FileSchema();
  FileSchema( const FileSchema & );
  FileSchema & operator = ( const FileSchema & );

  std::vector<const FilePart*> m_parts[ Element + 1 ];
};

//----------------------------------------------------------------------

struct FieldIO {
  const FieldBase * m_field ;
  unsigned          m_offset ;
  int               m_var_index ;
  const FilePart  * m_part ;
};


class FileOutput {
public:
  ~FileOutput();

  /** Create an output file for a collection of fields. */
  FileOutput( const FileSchema & ,
              const BulkData & ,
              const std::string & arg_file_path ,
              const std::string & arg_title ,
              const bool          arg_storage_double ,
              const std::vector< const FieldBase * > & ,
              const int * const arg_processor = NULL );

  /** Write a snapshot of field values */
  void write( double );

  const FileSchema & m_schema ;
  const BulkData       & m_mesh ;

  int exo_id() const { return m_exo_id ; }
  int exo_step() const { return m_counter ; }

  const std::vector<int> & global_counts() const ;
  const std::vector< FieldIO > & field_node_universal() const
    { return m_field_node_universal ; }

  const std::vector< FieldIO > & field_elem() const
    { return m_field_elem ; }

private:
  FileOutput();
  FileOutput( const FileOutput & );
  FileOutput & operator = ( const FileOutput & );

  int m_exo_id ;
  int m_counter ;
  int m_max_buffer ;

  std::vector< FieldIO > m_field_node_universal ;
  std::vector< FieldIO > m_field_elem ;
  std::vector<int> m_global_counts ;
};

//----------------------------------------------------------------------

class FileInput {
public:
  ~FileInput();

  FileInput( const FileSchema & , BulkData & ,
             const std::string & arg_file_path ,
             const std::vector< const FieldBase * > & );

  double read();

  const FileSchema & m_schema ;
        BulkData       & m_mesh ;

  int exo_id() const { return m_exo_id ; }
  int exo_step() const { return m_counter ; }

private:
  FileInput();
  FileInput( const FileInput & );
  FileInput & operator = ( const FileInput & );

  int m_exo_id ;
  int m_counter ;
  int m_max_buffer ;

  std::vector< FieldIO > m_field_node_universal ;
  std::vector< FieldIO > m_field_elem ;
};

} // namespace exodus
} // namespace phdmesh

//----------------------------------------------------------------------

#endif

