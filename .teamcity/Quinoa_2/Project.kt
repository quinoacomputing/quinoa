package Quinoa_2

import Quinoa_2.vcsRoots.*
import Quinoa_2.vcsRoots.Quinoa_2_GitGithubComQuinoacomputingQuinoaGitRefsHeadsMaster
import jetbrains.buildServer.configs.kotlin.v10.*
import jetbrains.buildServer.configs.kotlin.v10.Project
import jetbrains.buildServer.configs.kotlin.v10.projectFeatures.VersionedSettings
import jetbrains.buildServer.configs.kotlin.v10.projectFeatures.VersionedSettings.*
import jetbrains.buildServer.configs.kotlin.v10.projectFeatures.versionedSettings

object Project : Project({
    uuid = "e2dcb749-3b0c-44d8-ba6f-3e84b68b0fb0"
    extId = "Quinoa_2"
    parentId = "_Root"
    name = "Quinoa"

    vcsRoot(Quinoa_2_GitGithubComQuinoacomputingQuinoaGitRefsHeadsMaster)

    features {
        versionedSettings {
            id = "PROJECT_EXT_5"
            mode = VersionedSettings.Mode.ENABLED
            buildSettingsMode = VersionedSettings.BuildSettingsMode.PREFER_SETTINGS_FROM_VCS
            rootExtId = Quinoa_2_GitGithubComQuinoacomputingQuinoaGitRefsHeadsMaster.extId
            showChanges = false
            settingsFormat = VersionedSettings.Format.KOTLIN
        }
    }
})
