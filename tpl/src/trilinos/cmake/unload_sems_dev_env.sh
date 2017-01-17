#
# Unload the current Trilinos SEMS development environment
#
# USAGE:
#
#   $ source unload_sems_dev_env.sh
#
# If no Trilinos SEMS Dev Env is currently loaded, then this is a do nothing
# script.
#
# NOTE: If the modules don't unload correctly (e.g. because the
# sems-clang/<version> module loaded a GCC module), then consider calling:
#
#   $ module purge
#   $ unset TRILINOS_SEMS_DEV_ENV_LOADED
#
# to clean house (but this will also unload the sems-env module which will
# result in none of the other sems-xxx modules being visible).
#

# Other defaults (must be the same as in load_sems_dev_env.sh)
sems_python_and_version_default=sems-python/2.7.9
sems_boost_and_version_default=sems-boost/1.55.0/base

#
# A) Exit if no Trilinos SEMS Dev Env is currently loaded
#

if [ "$TRILINOS_SEMS_DEV_ENV_LOADED" == "" ] ; then
  if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
    echo "TRILINOS_SEMS_DEV_ENV_LOADED='', nothing to unload!"
  fi
  return 0
fi

#
# B) Get the loaded dev env info to unload the correct modules
#

idx=0
for str in $TRILINOS_SEMS_DEV_ENV_LOADED; do
  #echo $str
  sems_dev_env_array[$idx]=$str
  idx=$(expr $idx+1)
done
sems_compiler_and_version_unload=${sems_dev_env_array[0]}
sems_mpi_and_version_unload=${sems_dev_env_array[1]}
sems_cmake_and_version_unload=${sems_dev_env_array[2]}
#echo "sems_compiler_and_version_unload = $sems_compiler_and_version_unload"
#echo "sems_mpi_and_version_unload = $sems_mpi_and_version_unload"
if [ "$TRILINOS_SEMS_DEV_ENV_VERBOSE" == "1" ] ; then
  echo "Unloading TRILINOS_SEMS_DEV_ENV_LOADED =" \
    "'$sems_compiler_and_version_unload" \
    "$sems_mpi_and_version_unload" \
    "$sems_cmake_and_version_unload'!"
fi

#
# C) Unload the modules (in the reverse order)
#

module unload sems-superlu/4.3/base
module unload sems-scotch/6.0.3/parallel 
module unload sems-parmetis/4.0.3/parallel 
module unload sems-netcdf/4.3.2/parallel 
module unload sems-hdf5/1.8.12/parallel 
module unload sems-zlib/1.2.8/base 
module unload $sems_boost_and_version_default
module unload $sems_mpi_and_version_unload
module unload $sems_compiler_and_version_unload
module unload $sems_cmake_and_version_unload
module unload $sems_python_and_version_default

if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
  module list
fi

#
# D) Wipe out the rememberd Trilinos SEMS Dev Env
#

export TRILINOS_SEMS_DEV_ENV_LOADED=
