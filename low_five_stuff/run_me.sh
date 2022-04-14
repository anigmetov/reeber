#!/bin/sh

set -e

export HDF5_PLUGIN_PATH=/home/narn/code/LowFive/build_debug/src
export HDF5_VOL_CONNECTOR="lowfive under_vol=0;under_info={};"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_PLUGIN_PATH

make amr_connected_components_float > /dev/null 2>&1

rm -f log_*.txt
rm -f Backtrace.*

mpirun -n 1 -l ./amr_connected_components_float -i 0.5 -a -x 1 -f "/level_0/dark_matter_density" -n inputs_float.h5 /dev/null >qq 2>&1
