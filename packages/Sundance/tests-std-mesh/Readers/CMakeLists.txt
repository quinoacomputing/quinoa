# CMake tests specification 



IF (TPL_ENABLE_ExodusII)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
        partitionExo
        SOURCES partitionExo.cpp
        COMM serial mpi
)




TRIBITS_COPY_FILES_TO_BINARY_DIR(SundanceMeshCopyFiles
  DEST_FILES wheel.exo
  EXEDEPS partitionExo
  )

ENDIF()




