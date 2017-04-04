package Quinoa_2_Mac.buildTypes

import Quinoa_2_Mac.buildParams.*
import jetbrains.buildServer.configs.kotlin.v10.*

class Quinoa_2_Mac_Build( bp: BuildParams ) : BuildType({

    template(Quinoa_2_Mac.buildTypes.Quinoa_2_Mac_Matrix)

    uuid = "7df011b4-4795-4b89-9d03-aba1b1cb53f7_$paramToId"
    extId = "Quinoa_2_Mac_Build_$paramToId"
    name = "${bp.buildtype.toString()}, ${bp.compiler.toString()}"
    description = "Mac matrix build instance"

    params {
        param("buildtype", bp.buildtype.toString())
        param("compiler", bp.compiler.toString())
    }
})
