
TRIBITS_ADD_EXECUTABLE(
  AdvanceDependencyExample
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  AdvDepFileCopy
  SOURCE_FILES inputs.xml
  EXEDEPS AdvanceDependencyExample
)
