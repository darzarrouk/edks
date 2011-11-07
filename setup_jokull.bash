# initializing Pythia
source /home/geomod/dev/pythia/setup_jokull.sh

# INTEL COMPILER
#export LD_LIBRARY_PATH=/opt/intel/lib/intel64:${LD_LIBRARY_PATH}


# MPI_EDKS and EDKS python wrapper
export EDKS_HOME='/home/geomod/dev/edks'
export PYTHONPATH=${EDKS_HOME}/MPI_EDKS:$PYTHONPATH
export PATH=${EDKS_HOME}/MPI_EDKS:$PATH
export EDKS_BIN=${EDKS_HOME}/bin

# NUMBER OF CORES FOR OPENMPI
export OMP_NUM_THREADS=12

# MPI4PY
mpi4py_HOME=/home/geomod/dev/mpi4py-1.2.2/glnx64/lib64/python2.6/site-packages
export PYTHONPATH=${mpi4py_HOME}:${PYTHONPATH}

# ICM package
export PYTHONPATH=/home/geomod/dev/ICM/src:${PYTHONPATH}

