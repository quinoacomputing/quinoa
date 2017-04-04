package Quinoa_2_Mac.buildTypes

import jetbrains.buildServer.configs.kotlin.v10.*

object Quinoa_2_Mac_MacBuildFromTmp : BuildType({
    template(Quinoa_2_Mac.buildTypes.Quinoa_2_Mac_Matrix)
    uuid = "7df011b4-4795-4b89-9d03-aba1b1cb53f7"
    extId = "Quinoa_2_Mac_MacBuildFromTmp"
    name = "MacBuildFromTmp"

    params {
        param("buildtype", "blah")
        param("compiler", "blah")
    }
})
