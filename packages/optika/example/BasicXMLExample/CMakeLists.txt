
TRIBITS_ADD_EXECUTABLE(
  BasicXMLExample
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  BasicExampleFileCopy
  SOURCE_FILES inputs.xml
  EXEDEPS BasicXMLExample
)

