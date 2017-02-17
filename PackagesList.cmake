# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


#
# Define the Trilinos packages
#
TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  Gtest                 commonTools/gtest                 SS
  ThreadPool            packages/ThreadPool               PS # Depends on Pthreads
  Kokkos                packages/kokkos                   PS
  Teuchos               packages/teuchos                  PS
  RTOp                  packages/rtop                     PS
  Sacado                packages/sacado                   PS
  Epetra                packages/epetra                   PS
  SCORECpcu             SCOREC/pcu                        SS
  SCORECgmi             SCOREC/gmi                        SS
  SCORECgmi_sim         SCOREC/gmi_sim                    SS
  SCORECapf             SCOREC/apf                        SS
  SCORECapf_sim         SCOREC/apf_sim                    SS
  SCORECmds             SCOREC/mds                        SS
  SCORECparma           SCOREC/parma                      SS
  SCORECspr             SCOREC/spr                        SS
  Zoltan                packages/zoltan                   PS
  Shards                packages/shards                   PS
  GlobiPack             packages/globipack                PS
  Triutils              packages/triutils                 PS
  Tpetra                packages/tpetra                   PS
  EpetraExt             packages/epetraext                PS
  Domi                  packages/domi                     EX
  Thyra                 packages/thyra                    PS
  Xpetra                packages/xpetra                   PS
  OptiPack              packages/optipack                 PS
  Isorropia             packages/isorropia                PS
  Pliris                packages/pliris                   PS
  Claps                 packages/claps                    EX
  AztecOO               packages/aztecoo                  PS
  Galeri                packages/galeri                   PS
  Amesos                packages/amesos                   PS
  Pamgen                packages/pamgen                   PS
  Zoltan2               packages/zoltan2                  SS
  Ifpack                packages/ifpack                   PS
  ML                    packages/ml                       PS
  Belos                 packages/belos                    PS
  ShyLU                 packages/shylu                    SS
  Amesos2               packages/amesos2                  SS
  SEACAS                packages/seacas                   SS # Depends on netcdf, optionally hdf5, xdmf, pamgen
  Trios                 packages/trios                    EX #temporary
  Komplex               packages/komplex                  PS
  Anasazi               packages/anasazi                  PS
  Ifpack2               packages/ifpack2                  PS
  Stratimikos           packages/stratimikos              PS
  FEI                   packages/fei                      PS
  Teko                  packages/teko                     SS
  TriKota               packages/TriKota                  SS
  Intrepid              packages/intrepid                 PS
  Intrepid2             packages/intrepid2                SS
  STK                   packages/stk                      SS # Depends on boost
  SCORECapf_zoltan      SCOREC/zoltan                     SS
  SCORECapf_stk         SCOREC/stk                        SS
  SCORECma              SCOREC/ma                         SS
  SCORECpumi            SCOREC/pumi                       SS
  SCOREC                SCOREC                            SS
  Phalanx               packages/phalanx                  SS
  NOX                   packages/nox                      PS
  Moertel               packages/moertel                  PS
  MueLu                 packages/muelu                    SS
  Rythmos               packages/rythmos                  PS
  Tempus                tempus                            ST
  MOOCHO                packages/moocho                   ST
  Stokhos               packages/stokhos                  SS
  ROL                   packages/rol                      SS
  Piro                  packages/piro                     SS
  Panzer                packages/panzer                   SS
  Sundance              packages/Sundance                 SS # Could be PS based on deps (BUG: 4669)
  CTrilinos             packages/CTrilinos                SS # Switched to SS to speed up checkin testing
  ForTrilinos           packages/ForTrilinos              EX
  PyTrilinos            packages/PyTrilinos               SS
  WebTrilinos           packages/WebTrilinos              EX # Should be SS
  NewPackage            packages/new_package              EX # Should be SS
  Optika		packages/optika		          EX
  Mesquite              packages/mesquite                 PS
  MeshingGenie          packages/meshinggenie             EX
  TrilinosCouplings     packages/trilinoscouplings        SS
  Pike                  packages/pike                     SS
  xSDKTrilinos          packages/xSDKTrilinos             SS
  )

# Allow builds even if some packages are missing

TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCOREC)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECgmi)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECgmi_sim)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECpcu)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_sim)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECmds)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECparma)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECspr)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_stk)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_zoltan)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECma)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECpumi)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Tempus)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(MOOCHO)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Sundance)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(CTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(ForTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Optika)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Mesquite)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(WebTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(xSDKTrilinos)

#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(MOOCHO Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Phalanx Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(PyTrilinos Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Sundance Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Tpetra Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Ifpack2 Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(TriKota Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Pamgen Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(STK Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(SEACAS Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Anasazi Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Zoltan Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Isorropia Windows)

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Teko Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Mesquite AIX)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Trios Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Panzer Windows)
