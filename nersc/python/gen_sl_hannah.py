#!/usr/bin/env python3

# /project/projectdirs/m2848/hannah_Nyx/nbody/1species_vog_256_refined2level_8crit/plotfiles

hannah_plt_files = [ "/project/projectdirs/m2848/hannah_Nyx/nbody/1species_vog_256_refined2level_8crit/plotfiles/plt00719" ]
# hannah_plt_files = [ "/project/projectdirs/m2848/hannah_Nyx/nbody/1species_vog_256_refined2level_8crit/plotfiles/plt00488" ]

hannah_plt_nums = [ "00719" ]
# hannah_plt_nums = [ "00488" ]


def write_sl(file_name, plt_file_name, plt_num, n_cores, max_minutes, n_run):
    mpi_per_node = 68
    n_nodes = 1 + n_cores // mpi_per_node
    c_knl = 272 // mpi_per_node
    rho = "81.66"
    theta = "10"

    npy_fname = "plt{0}".format(plt_num)
    # dgm_fname = "$SCRATCH/output/my_answers_knl/amr_si_{0}_blocks_{1}_rho_{2}_diagram.txt".format(npy_fname, n_cores, rho)
    # intl_fname = "$SCRATCH/output/my_answers_knl/amr_si_{0}_blocks_{1}_rho_{2}_theta_{3}_integral-non-absolute.txt".format(npy_fname, n_cores, rho, theta)
    # mt_fname = "$SCRATCH/output/my_answers_knl/amr_si_plt_num_{0}_blocks_{1}.mt".format(plt_num, n_cores)
    dgm_fname = "none"
    intl_fname = "none"
    mt_fname = "none"
    log_fname = "log_cc_cori_knl_hannah_plt_num_{0}_cores_{1}_run_{2}_rho_{3}_amr.txt".format(plt_num, n_cores, n_run, rho)

    script_template = """#!/bin/bash -l

#SBATCH -C knl
#SBATCH -A m2848
#SBATCH -q premium
#SBATCH -N {0}
#SBATCH -t 00:{1}:00
#SBATCH -L SCRATCH
#SBATCH -o {2}
#SBATCH -e {2}

srun -n {3} -c {10}  --cpu-bind=cores ./amr_connected_components_float -f particle_mass_density -r 10 --plotfile -i {4} -n -w -x {5} {6} {7} {8} {9}\n"""

    with open(file_name, 'w') as sl_file:
        sl_file.write(script_template.format(n_nodes, max_minutes, log_fname, n_cores, rho, theta, plt_file_name, mt_fname, dgm_fname, intl_fname, c_knl))


if __name__ == "__main__":
    sbatch_file_name = "queue_all_hannah.sh"
    max_minutes = 2
    with open(sbatch_file_name, 'w') as sbatch_file:
        sbatch_file.write('#!/bin/sh\n\n')
        for plt_file_name, plt_num in zip(hannah_plt_files, hannah_plt_nums):
            # for n_run in range(0,3):
            for n_run in [50]:
                for n_cores in [4096,8192]:
                # for n_cores in [2048, 4096, 4096 * 2, 4096 * 4, 3045, 6089, 6089*2, 6089 * 4]:
                # for n_cores in [761, 1522, 3045, 6089, 6089 * 2]:
                # for n_cores in [761]:
                # for n_cores in [6089*2, 6089]:
                    file_name = "run_cc_hannah_{0}_{1}_{2}_amr.sl".format(plt_num, n_cores, n_run)
                    write_sl(file_name, plt_file_name, plt_num, n_cores, max_minutes, n_run)
                    sbatch_file.write("sbatch {0}\n".format(file_name))
