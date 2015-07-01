
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(DualScopeSet)

#
# A) Define the package
#

TRIBITS_PACKAGE(WebTrilinos  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(c++)


#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#TRIBITS_EXCLUDE_FILES(
#  util
#  etc
#  test/tmp
#  doc/UsersGuide
#  doc/DevelopersGuide
#  doc/MLAPI
#  python
#  doc/DoxyfileWeb
#  doc/build_docs
#  doc/ml-logo.eps
#  doc/ml-logo.jpg
#  doc/sc2000.ps.gz
#  examples/Makefile-common.include
#  examples/Maxwell/ml_periodic_max.c
#  examples/Other/ml_read_complexmaxwell.c
#  examples/Other/ml_read_maxwell.c
#  examples/Other/ml_star2d.c
#  examples/Other/new_readex.c
#  examples/Other/oldml_readex.c
#  examples/Other/seg_readex.c
#  examples/README.AddingExamples
#  examples/RefMaxwell
#  examples/RefMaxwell/rpc.cpp
#  src/Coarsen/README
#  src/Main/ml_v_cycle.c
#  src/Smoother/README
#  src/Utils/jmpilib.c
#  src/Utils/jostle.h
#  src/Utils/ml_vampir.c
#  src/Utils/tumi.c
#  test/README.runtests
#  test/Zoltan/cxx_main_simple.cpp
#  test/scripts
#  test/scripts/run-tests.sh
#  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
