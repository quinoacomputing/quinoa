# Generate Mesquite_all_headers.hpp from list of header files

separate_arguments( headers )
file( WRITE ${output_file} "#ifndef MESQUITE_ALL_HEADERS_HPP\n#define MESQUITE_ALL_HEADERS_HPP\n" )
foreach( header ${headers} )
  string( REGEX MATCH "[^/]+\\.h(pp)?" file ${header} )
  file( APPEND ${output_file} "#include \"${file}\"\n")
endforeach()
file( APPEND ${output_file} "#endif\n" )
