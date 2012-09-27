# initializing Pythia
source /home/ortega/tools/pythia0.8/setup.sh

# MPI_EDKS and EDKS python wrapper
export EDKS_HOME='/home/ortega/dev/edks'
export PYTHONPATH=${EDKS_HOME}/MPI_EDKS:$PYTHONPATH
export PATH=${EDKS_HOME}/MPI_EDKS:$PATH
export EDKS_BIN=${EDKS_HOME}/bin

# NUMBER OF CORES FOR OPENMPI
export OMP_NUM_THREADS=11

# ICM package
export PYTHONPATH=/home/ortega/dev/ICM/src:${PYTHONPATH}

# fix bug that does not assign the current folder to the pythonpath
export PYTHONPATH=.:${PYTHONPATH}
