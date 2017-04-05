package Quinoa_2_Docker.buildTypes

import jetbrains.buildServer.configs.kotlin.v10.*

object Quinoa_2_Docker_Debian : BuildType({
    template(Quinoa_2_Docker.buildTypes.Quinoa_2_Docker_Image)
    uuid = "95961fba-6827-45d0-af17-46b9368f76d3"
    extId = "Quinoa_2_Docker_Debian"
    name = "Debian"
    description = "quinoacomputing/quinoa:debian@hub.docker.com"

    params {
        param("dockerfile", "Dockerfile.quinoa-build-debian")
        param("organization", "quinoacomputing")
        param("repository", "quinoa")
        param("tag", "debian")
        param("workdir", "docker")
    }
})
