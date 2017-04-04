package Quinoa_2_Linux.buildTypes

import jetbrains.buildServer.configs.kotlin.v10.*

object Quinoa_2_Linux_BuildFromMatrix : BuildType({
    template(Quinoa_2_Linux.buildTypes.Quinoa_2_Linux_Matrix)
    uuid = "1308360c-2059-48e8-8188-c06d1f15ecfb"
    extId = "Quinoa_2_Linux_BuildFromMatrix"
    name = "BuildFromMatrix"

    params {
        param("buildtype", "blah")
        param("compiler", "blah")
        param("mathlib", "blah")
        param("rngsse2", "blah")
        param("stdlibcpp", "blah")
        param("testu01", "blah")
    }
})
