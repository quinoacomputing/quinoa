INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE(
  ValidatorExample
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  ValidatorExampleFileCopy
  SOURCE_FILES inputs.xml
  EXEDEPS ValidatorExample
)

