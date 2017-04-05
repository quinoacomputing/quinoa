package Quinoa_2_Docker.buildTypes

import jetbrains.buildServer.configs.kotlin.v10.*

object Quinoa_2_Docker_Alpine : BuildType({
    template(Quinoa_2_Docker.buildTypes.Quinoa_2_Docker_Image)
    uuid = "da2c6f8c-ebc4-4a09-aa57-8a004eca318f"
    extId = "Quinoa_2_Docker_Alpine"
    name = "Alpine"
    description = "quinoacomputing/quinoa:alpine@hub.docker.com"

    params {
        param("dockerfile", "Dockerfile.quinoa-build-alpine")
        param("organization", "quinoacomputing")
        param("repository", "quinoa")
        param("tag", "alpine")
        param("workdir", "docker")
    }
})
