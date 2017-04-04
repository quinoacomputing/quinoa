package Quinoa_2_Mac

import Quinoa_2_Mac.buildTypes.*
import jetbrains.buildServer.configs.kotlin.v10.*
import jetbrains.buildServer.configs.kotlin.v10.Project

object Project : Project({
    uuid = "7b46bd08-3540-4fa1-ab86-d238159396f3"
    extId = "Quinoa_2_Mac"
    parentId = "Quinoa_2"
    name = "Mac"
    description = "Mac builds"

    buildType(Quinoa_2_Mac_MacBuildFromTmp)

    template(Quinoa_2_Mac_Matrix)
})
