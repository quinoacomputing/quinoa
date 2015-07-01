
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

# dummy change to test secondary git repo
# another dummy change

#
# A) Define the package
#

#TRIBITS_PACKAGE(Playa  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#

# ------- nothing to do here

#
# C) Add the libraries, tests, and examples
#

MESSAGE("adding Playa/src")
ADD_SUBDIRECTORY(src)
MESSAGE("adding Playa/tests")
TRIBITS_ADD_TEST_DIRECTORIES(tests)
MESSAGE("adding Playa/examples")
TRIBITS_ADD_TEST_DIRECTORIES(examples)


#
# Exclude files for source package.
#

#TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

#TRIBITS_PACKAGE_POSTPROCESS()
