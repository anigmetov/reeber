#!/usr/bin/env python3

def write_sl(file_name, plt_num, n_cores, max_minutes):
    mpi_per_node =  68
    n_nodes = 1 + n_cores // (mpi_per_node)
    queue_name = "regular"
    c_knl = (4 * 68) // mpi_per_node
    rho = "81.66"
    min_cells = "10"

    npy_fname = "plt{0}".format(plt_num)
    dgm_fname = "none"
    intl_fname = "$SCRATCH/output/my_answers_knl/zarija/new/{0}/integral_{0}_blocks_{1}_rho_{2}".format(npy_fname, n_cores, rho)
    npy_full_fname = "Pltfiles/{0}".format(npy_fname)
    mt_fname = "none"
    log_fname = "log_cc_print_halo_plt_num_{0}_cores_{1}_amr.txt".format(plt_num, n_cores)

    script_template = """#!/bin/bash -l

#SBATCH -C knl
#SBATCH -A m2848
#SBATCH -q {11}
#SBATCH -N {0}
#SBATCH -t 00:{1}:00
#SBATCH -L SCRATCH
#SBATCH -o {2}
#SBATCH -e {2}

srun -n {3} -c {10}  --cpu-bind=cores ./amr_connected_components_float -f density,particle_mass_density,xmom,ymom,zmom --plotfile -i {4} -n -w -x {5} {6} {7} {8} {9}\n"""

    with open(file_name, 'w') as sl_file:
        sl_file.write(script_template.format(n_nodes, max_minutes, log_fname, n_cores, rho, min_cells, 
            npy_full_fname, mt_fname, dgm_fname, intl_fname, c_knl, queue_name))

# IC     TREECOOL_ONO_INST_HeIIHomo  knl.sl    plt00475  plt00552  plt00632  plt00712  plt00788  plt00881  plt00992  plt01121  plt01272  probin  zhi_2048
# Slurm  inputs                      plt00300  plt00509  plt00594  plt00678  plt00749  plt00831  plt00934  plt01055  plt01198  plt01354  runlog



if __name__ == "__main__":
    sbatch_file_name = "queue_all_print_halo.sh"
    max_minutes = 14

    plt_nums = [ "00300", "00475", "00509", "00552", "00594", "00632", "00678",
                 "00712", "00749", "00788", "00831", "00881", "00934", "00992",
                 "01055", "01121", "01198", "01272", "01354" ]

    plt_nums = [ "00632", "00678",
                 "00712", "00749", "00788", "00831", "00881", "00934",
                 "01055", "01121", "01198", "01272", "01354" ]

 
    with open(sbatch_file_name, 'w') as sbatch_file:
        sbatch_file.write('#!/bin/sh\n\n')
        for plt_num in plt_nums:
        # for plt_num in ["01198", "01272"]:
        # for plt_num in ["01121", "01055", "00992", "00881", "00831", "00788", "00749", "00712"]:
            for n_cores in [32 * 2048]:
                file_name = "run_cc_big_print_halo_{0}_{1}_amr.sl".format(plt_num, n_cores)
                write_sl(file_name, plt_num, n_cores, max_minutes)
                sbatch_file.write("mkdir -p new_answers/plt{0}\n".format(plt_num))
                sbatch_file.write("sbatch {0}\n".format(file_name))

