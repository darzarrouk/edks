# initializing Pythia
source /Users/ortega/tools/pythia-0.8/setup.sh

# INTEL COMPILER
#export LD_LIBRARY_PATH=/opt/intel/lib/intel64:${LD_LIBRARY_PATH}


# MPI_EDKS and EDKS python wrapper
export EDKS_HOME=${HOME}/dev/edks
export PYTHONPATH=${EDKS_HOME}/MPI_EDKS:$PYTHONPATH
export PATH=${EDKS_HOME}/MPI_EDKS:$PATH
export EDKS_BIN=${EDKS_HOME}/bin

# MPI4PY
mpi4py_HOME=/Users/ortega/tools/mpi4py-1.2.2/darwin/lib/python
export PYTHONPATH=${mpi4py_HOME}:${PYTHONPATH}

# ICM package
export PYTHONPATH=${HOME}/dev/ICM/src:${PYTHONPATH}

