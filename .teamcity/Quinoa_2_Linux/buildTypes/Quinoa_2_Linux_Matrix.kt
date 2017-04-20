package Quinoa_2_Linux.buildTypes

import jetbrains.buildServer.configs.kotlin.v10.*
import jetbrains.buildServer.configs.kotlin.v10.buildSteps.ScriptBuildStep
import jetbrains.buildServer.configs.kotlin.v10.buildSteps.ScriptBuildStep.*
import jetbrains.buildServer.configs.kotlin.v10.buildSteps.script
import jetbrains.buildServer.configs.kotlin.v10.triggers.VcsTrigger
import jetbrains.buildServer.configs.kotlin.v10.triggers.VcsTrigger.*
import jetbrains.buildServer.configs.kotlin.v10.triggers.vcs

object Quinoa_2_Linux_Matrix : Template({
    uuid = "dfd6c8ef-72e2-4e44-be66-fea7a62d455e"
    extId = "Quinoa_2_Linux_Matrix"
    name = "Matrix"

    vcs {
        root(Quinoa_2.vcsRoots.Quinoa_2_GitGithubComQuinoacomputingQuinoaGitRefsHeadsMaster)

    }

    val stepPrefix = """
      . ${'$'}SPACK_ROOT/share/spack/setup-env.sh
      [ %compiler% == clang ] && module load clang/latest openmpi/1.10.2/clang/latest hdf5-1.10.0-patch1-gcc-4.8.5-mmtlfty netcdf-4.4.1-gcc-4.8.5-5xen4a5 hypre-2.10.1-gcc-4.8.5-beuxbxv
      [ %compiler% == gnu ] && module load openmpi-2.0.1-gcc-4.8.5-jv7w2de hdf5-1.10.0-patch1-gcc-4.8.5-mmtlfty netcdf-4.4.1-gcc-4.8.5-5xen4a5 hypre-2.10.1-gcc-4.8.5-beuxbxv
      [ %compiler% == intel ] && module load intel/2018 mpi hdf5/intel netcdf/intel hypre/intel
      [ %mathlib% == mkl ] && module load mkl/2018
      [[ %mathlib% == lapack && %compiler% == intel ]] && module load lapack/intel
      [[ %mathlib% == lapack && %compiler% != intel ]] && module load netlib-lapack-3.6.1-gcc-4.8.5-snwxnfw
      [ %rngsse2% == true ] && module load rngsse2
      [ %testu01% == true ] && module load testu01
      module load charm/%compiler%-%stdlibcpp% h5part/%compiler% trilinos/%compiler%-%stdlibcpp%/%mathlib%
      module load pugixml pegtl pstreams boost-1.61.0-gcc-4.8.5-q2hywin gmsh-2.12.0-gcc-4.8.5-p3vpjfb random123 tut cartesian_product numdiff libc++
      module list""".trimIndent()

    steps {
        script {
            name = "Verify commit"
            id = "RUNNER_17"
            scriptContent = """/ccs/opt/git/bin/git verify-commit %build.vcs.number% 2>&1 | grep "Good signature""""
        }
        script {
            name = "Build code"
            id = "RUNNER_18"
            scriptContent = """
                ${stepPrefix}
                rm -rf build && mkdir build && cd build
                cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=%buildtype% -DSTDLIBCPP=%stdlibcpp% -DCMAKE_DISABLE_FIND_PACKAGE_RNGSSE2=!%rngsse2% -DCMAKE_DISABLE_FIND_PACKAGE_TestU01=!%testu01% -DCMAKE_CXX_FLAGS=-Werror ../src
                make -j16
            """.trimIndent()
        }
        script {
            name = "Run tests"
            id = "RUNNER_19"
            workingDir = "build"
            scriptContent = """
                ${stepPrefix}
                ../script/run_tests.sh 16
            """.trimIndent()
        }
    }

    requirements {
        equals("teamcity.agent.jvm.os.name", "Linux", "RQ_19")
        contains("teamcity.agent.name", "ccscs", "RQ_20")
    }

    triggers {
        vcs {
            id = "vcsTrigger"
            triggerRules = """
                +:.
                -:comment=\[ci skip\]:**
                -:comment=\[skip ci\]:**
            """.trimIndent()
            branchFilter = "+:<default>"
            perCheckinTriggering = true
            groupCheckinsByCommitter = true
        }
    }
})
