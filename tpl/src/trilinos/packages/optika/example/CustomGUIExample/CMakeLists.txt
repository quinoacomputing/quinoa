INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE(
  CustomGUIExample
  SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR( CustomGUIFileCopies
  DEST_FILES
  myStyleSheet.qss myIcon.png
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/example/CustomGUIExample
  EXEDEPS CustomGUIExample
  )
