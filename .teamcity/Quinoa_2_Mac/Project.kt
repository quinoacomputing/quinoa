package Quinoa_2_Mac

import Quinoa_2_Mac.buildTypes.*
import Quinoa_2_Mac.buildParams.*
import jetbrains.buildServer.configs.kotlin.v10.*
import jetbrains.buildServer.configs.kotlin.v10.Project

object Project : Project({
    uuid = "7b46bd08-3540-4fa1-ab86-d238159396f3"
    extId = "Quinoa_2_Mac"
    parentId = "Quinoa_2"
    name = "Mac"
    description = "Mac builds"

    template(Quinoa_2_Mac_Matrix)
    val allBuilds = mutableListOf< BuildParams >()

    // Generate matrix with all possible combinations of build parameters,
    // defined in package buildParams.
    CmakeBuildType.values().forEach{ b ->
      Compiler.values().forEach{ c ->
         allBuilds.add( BuildParams(b,c) )
      }
    }

    val builds = mutableListOf< BuildParams >()

    // Optionally exclude some builds
    allBuilds.forEach{ b ->
      builds.add( b );
    }

    // Generate TeamCity builds
    builds.forEach{ buildType( Quinoa_2_Mac_Build(it) ) }
})
