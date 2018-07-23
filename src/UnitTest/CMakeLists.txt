cmake_minimum_required(VERSION 2.8.5)

project(UnitTest CXX)

# Include function for adding Charm++ modules
include(charm)

if (HAS_MKL)
  set(TestMKLBetaMethod "tests/Control/Options/TestMKLBetaMethod.C")
  set(TestMKLGammaMethod "tests/Control/Options/TestMKLGammaMethod.C")
  set(TestMKLGaussianMethod "tests/Control/Options/TestMKLGaussianMethod.C")
  set(TestMKLUniformMethod "tests/Control/Options/TestMKLUniformMethod.C")
  set(TestMKLRNG "tests/RNG/TestMKLRNG.C")
endif()

if(HAS_RNGSSE2)
  set(TestRNGSSE "tests/RNG/TestRNGSSE.C")
endif()

add_library(UnitTest
            Assessment.C
            TUTSuite.C
            TUTTest.C
            tests/Base/TestContainerUtil.C
            tests/Base/TestData.C
            tests/Base/TestException.C
            tests/Base/TestExceptionMPI.C
            tests/Base/TestFactory.C
            tests/Base/TestFlip_map.C
            tests/Base/TestHas.C
            tests/Base/TestPrint.C
            tests/Base/TestProcessControl.C
            tests/Base/TestPUPUtil.C
            tests/Base/TestReader.C
            tests/Base/TestStrConvUtil.C
            tests/Base/TestTaggedTuple.C
            tests/Base/TestTimer.C
            tests/Base/TestVector.C
            tests/Base/TestWriter.C
            ${TestMKLUniformMethod}
            ${TestMKLGaussianMethod}
            ${TestMKLBetaMethod}
            ${TestMKLGammaMethod}
            tests/Control/Options/TestRNG.C
            tests/Control/TestControl.C
            tests/Control/TestFileParser.C
            tests/Control/TestStringParser.C
            tests/Control/TestSystemComponents.C
            tests/Control/TestToggle.C
            tests/Inciter/TestScheme.C
            tests/Inciter/AMR/TestError.C
            tests/IO/TestExodusIIMeshReader.C
            tests/IO/TestMesh.C
            tests/IO/TestMeshReader.C
            tests/LoadBalance/TestLinearMap.C
            tests/LoadBalance/TestLoadDistributor.C
            tests/LoadBalance/TestUnsMeshMap.C
            tests/Mesh/TestAround.C
            tests/Mesh/TestDerivedData.C
            tests/Mesh/TestGradients.C
            tests/Mesh/TestReorder.C
            ${TestMKLRNG}
            ${TestRNGSSE}
            tests/RNG/TestRNG.C
            tests/RNG/TestRandom123.C
)

addCharmModule( "tuttest" "UnitTest" )
addCharmModule( "tutsuite" "UnitTest" )
addCharmModule( "charmchild" "UnitTest" )
addCharmModule( "charmtimer" "UnitTest" )
addCharmModule( "testarray" "UnitTest" )
addCharmModule( "migrated_base" "UnitTest" )
addCharmModule( "migrated_inciter" "UnitTest" )

# Add extra dependency of UnitTest on unittestCharmModule. This is required as a
# dependency of UnitTest, TUTSuite, refers to the main Charm++ proxy defined in
# the Charm++ module unittest (in Main/UnitTest.C).
add_dependencies("UnitTest" "unittestCharmModule")

# Add extra dependencies of UnitTest on those CharmModules that are required for
# testing Inciter.
add_dependencies("UnitTest" "matcgCharmModule")
add_dependencies("UnitTest" "diagcgCharmModule")
add_dependencies("UnitTest" "distfctCharmModule")
add_dependencies("UnitTest" "dgCharmModule")
add_dependencies("UnitTest" "discretizationCharmModule")
add_dependencies("UnitTest" "solverCharmModule")
add_dependencies("UnitTest" "transporterCharmModule")

# Add extra dependency of TUTSuite on charm modules testing advanced array
# element placement using maps. This is required since TUTSuite spawns chare
# arrays testing advanced array element placement maps.
add_dependencies( "tutsuiteCharmModule" "linearmapCharmModule" )
add_dependencies( "tutsuiteCharmModule" "unsmeshmapCharmModule" )

set_target_properties(UnitTest PROPERTIES LIBRARY_OUTPUT_NAME quinoa_unittest)

INSTALL(TARGETS UnitTest
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
