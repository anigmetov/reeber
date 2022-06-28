#!/usr/bin/env python3

fred21_plt_files = [ "/global/cscratch1/sd/zarija/21Fred/L20_N1024_z3_plt00849",
                "/global/cscratch1/sd/zarija/21Fred/L20_N2048_z3_plt01846",
                "/global/cscratch1/sd/zarija/21Fred/L20_N512_z3_plt00414",
                "/global/cscratch1/sd/zarija/21Fred/L20_N256_z3_plt00253"]

# fred21_plt_files = [ "/global/cscratch1/sd/zarija/21Fred/L20_N2048_z3_plt01846" ]

fred21_plt_nums = [ 1024, 2048, 512, 256 ]
# fred21_plt_nums = [ 2048 ]


def write_sl(file_name, plt_file_name, plt_num, n_cores, max_minutes):
    mpi_per_node = 4*64
    n_nodes = 1 + n_cores // 68
    c_knl = 272 // mpi_per_node
    rho = "81.66"
    theta = "10"

    npy_fname = "plt{0}".format(plt_num)
    dgm_fname = "none"
    intl_fname = "my_answers_fred/intl_fred_zarija_num_{0}_cores_{1}".format(plt_num, n_cores)
    mt_fname = "none"
    log_fname = "log_cc_cori_knl_fred_plt_num_{0}_cores_{1}_amr.txt".format(plt_num, n_cores)

    script_template = """#!/bin/bash -l

#SBATCH -C knl
#SBATCH -q debug
#SBATCH -N {0}
#SBATCH -t 00:{1}:00
#SBATCH -L SCRATCH
#SBATCH -o {2}
#SBATCH -e {2}

srun -n {3} -c {10}  --cpu-bind=cores ./amr_connected_components_float --plotfile -i {4} -f particle_mass_density,density,xmom,ymom,zmom -n -w -x {5} {6} {7} {8} {9}\n"""

    with open(file_name, 'w') as sl_file:
        sl_file.write(script_template.format(n_nodes, max_minutes, log_fname, n_cores, rho, theta, plt_file_name, mt_fname, dgm_fname, intl_fname, c_knl))


if __name__ == "__main__":
    sbatch_file_name = "queue_all_fred.sh"
    max_minutes = 3
    with open(sbatch_file_name, 'w') as sbatch_file:
        sbatch_file.write('#!/bin/sh\n\n')
        for plt_file_name, plt_num in zip(fred21_plt_files, fred21_plt_nums):
            for n_cores in [1024]:
                file_name = "run_cc_fred_{0}_{1}_amr.sl".format(plt_num, n_cores)
                write_sl(file_name, plt_file_name, plt_num, n_cores, max_minutes)
                sbatch_file.write("sbatch {0}\n".format(file_name))
