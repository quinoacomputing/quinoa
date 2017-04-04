package Quinoa_2_Linux

import Quinoa_2_Linux.buildTypes.*
import Quinoa_2_Linux.buildParams.*
import jetbrains.buildServer.configs.kotlin.v10.*
import jetbrains.buildServer.configs.kotlin.v10.Project

object Project : Project({
    uuid = "19ac1646-b1ed-4d1f-bed9-a16af31a8062"
    extId = "Quinoa_2_Linux"
    parentId = "Quinoa_2"
    name = "Linux"
    description = "Linux builds"

    template(Quinoa_2_Linux_Matrix)

    val allBuilds = mutableListOf< BuildParams >()

    // Generate matrix with all possible combinations of build parameters,
    // defined in package buildParams.
    Compiler.values().forEach{ c ->
      StdLibC.values().forEach{ l ->
        MathLib.values().forEach{ m ->
          CmakeBuildType.values().forEach{ b ->
            for( r in listOf( true, false ) ) {
              for( t in listOf( true, false ) ) {
                allBuilds.add( BuildParams(b,c,m,l,r,t) )
              }
            }
          }
        }
      }
    }

    val builds = mutableListOf< BuildParams >()

    // Exclude some builds
    allBuilds.forEach{ b ->
      if ( !(b.compiler == Compiler.gnu && b.stdlibc == StdLibC.libc)) {
        builds.add( b );
      }
    }

    // Generate TeamCity builds
    builds.forEach{ buildType( QuinoaKotlin_Linux_Build(it) ) }
})
