package Quinoa_2.vcsRoots

import jetbrains.buildServer.configs.kotlin.v10.*
import jetbrains.buildServer.configs.kotlin.v10.vcs.GitVcsRoot

object Quinoa_2_GitGithubComQuinoacomputingQuinoaGitRefsHeadsMaster : GitVcsRoot({
    uuid = "fe8fc2bd-e317-494e-8018-baf9242fcf87"
    extId = "Quinoa_2_GitGithubComQuinoacomputingQuinoaGitRefsHeadsMaster"
    name = "git@github.com:quinoacomputing/quinoa.git#refs/heads/master"
    url = "git@github.com:quinoacomputing/quinoa.git"
    branchSpec = "+:refs/heads/(*)"
    checkoutSubmodules = GitVcsRoot.CheckoutSubmodules.IGNORE
    authMethod = defaultPrivateKey {
        userName = "git"
    }
})
