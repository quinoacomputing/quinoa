package Quinoa_2_Linux.buildTypes

import Quinoa_2_Linux.buildParams.*
import jetbrains.buildServer.configs.kotlin.v10.*

class Quinoa_2_Linux_Build( bp: BuildParams ) : BuildType({

    template(Quinoa_2_Linux.buildTypes.Quinoa_2_Linux_Matrix)

    val paramToId = bp.buildtype.toString().toExtId() +
                    bp.compiler.toString().toExtId() +
                    bp.mathlib.toString().toExtId() +
                    bp.stdlibc.toString().toExtId() +
                    bp.rngsse2.toString().toExtId() +
                    bp.testu01.toString().toExtId()

    uuid = "1308360c-2059-48e8-8188-c06d1f15ecfb_$paramToId"
    extId = "Quinoa_2_Linux_Build_$paramToId"
    name = "${bp.buildtype.toString()}, ${bp.compiler.toString()}, ${bp.mathlib.toString()}, ${bp.stdlibc.toString()}++, rngsse2=${bp.rngsse2}, testu01=${bp.testu01}"
    description = "Linux matrix build instance"

    params {
        param("buildtype", bp.buildtype.toString())
        param("compiler", bp.compiler.toString())
        param("mathlib", bp.mathlib.toString())
        param("stdlibcpp", bp.stdlibc.toString()+"++")
        param("rngsse2", bp.rngsse2.toString())
        param("testu01", bp.testu01.toString())
    }
})
