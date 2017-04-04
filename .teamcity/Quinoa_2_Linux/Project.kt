package Quinoa_2_Linux

import Quinoa_2_Linux.buildTypes.*
import jetbrains.buildServer.configs.kotlin.v10.*
import jetbrains.buildServer.configs.kotlin.v10.Project

object Project : Project({
    uuid = "19ac1646-b1ed-4d1f-bed9-a16af31a8062"
    extId = "Quinoa_2_Linux"
    parentId = "Quinoa_2"
    name = "Linux"
    description = "Linux builds"

    buildType(Quinoa_2_Linux_BuildFromMatrix)

    template(Quinoa_2_Linux_Matrix)
})
