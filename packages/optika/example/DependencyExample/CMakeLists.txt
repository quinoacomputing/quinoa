INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE(
  DependencyExample
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  DependencyExampleXMLCopy
  SOURCE_FILES inputs.xml
  EXEDEPS DependencyExample
)
