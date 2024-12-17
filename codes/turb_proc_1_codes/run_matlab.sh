#!/bin/bash
#SBATCH --job-name=isoplot
#SBATCH --time=2:00:00
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --account=geyink1
### Using more tasks because default memory is ~5GB per core
### 'shared' will share the node with other users
### 'parallel' use entire node (24,28,48, depends on node type)
### Try the script with --ntasks-per-node=1 and see what happens
#
#---------------------------------------------------------------------
# SLURM job script to run serial MATLAB
#---------------------------------------------------------------------
 
#ml matlab
#ml java
ml # confirm modules used
ml matlab
matlab -nodisplay -nosplash -nodesktop -r "calc_velfields_1"

#calc_velgrad_x;clear all;
#calc_velgrad_y;clear all;
#calc_velgrad_z;clear all;
#calc_vort;clear all;
#calc_conv_visc;clear all;
#calc_profile;"
echo "matlab exit code: $?"
