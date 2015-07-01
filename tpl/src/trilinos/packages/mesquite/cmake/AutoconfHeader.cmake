macro ( autoconf_header INFILE OUTFILE )

file( READ ${INFILE} autoconf_HEADER_IN )

string( REGEX MATCHALL "#undef ([^\n]*)" autoconf_VARS "${autoconf_HEADER_IN}" )
set( autoconf_HEADER "${autoconf_HEADER_IN}" )
foreach ( VAR ${autoconf_VARS} )
  string ( REGEX REPLACE "#undef (.*)" "\\1" VAR "${VAR}" )
  # CMake variables should usually be left undefined if they are blank or 0...
  if ( ${VAR} )
    string( REGEX REPLACE "#undef ${VAR}\n" "#define ${VAR} ${${VAR}}\n" autoconf_HEADER "${autoconf_HEADER}" )
  endif ( ${VAR} )
  # ... but always define version numbers, even if they have "0" as a value
#  if ( \"${VAR}\" MATCHES ".*VERSION.*" )
#    string( REGEX REPLACE "#undef ${VAR}\n" "#define ${VAR} ${${VAR}}\n" autoconf_HEADER "${autoconf_HEADER}" )
#  endif ( \"${VAR}\" MATCHES ".*VERSION.*" )
endforeach ( VAR )
string( CONFIGURE "${autoconf_HEADER}"  autoconf_HEADER_OUT )

if ( EXISTS "${OUTFILE}" )
  file( READ "${OUTFILE}" __autoconf_HEADER_PREV )
  if ( NOT "${autoconf_HEADER_OUT}" STREQUAL "${__autoconf_HEADER_PREV}" )
    file( WRITE "${OUTFILE}" "${autoconf_HEADER_OUT}" )
  endif ( NOT "${autoconf_HEADER_OUT}" STREQUAL "${__autoconf_HEADER_PREV}" )
else ( EXISTS "${OUTFILE}" )
  file( WRITE "${OUTFILE}" "${autoconf_HEADER_OUT}" )
endif ( EXISTS "${OUTFILE}" )

endmacro( autoconf_header )
