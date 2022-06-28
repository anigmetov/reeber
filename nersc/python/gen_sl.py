#!/usr/bin/env python3

def write_sl(file_name, plt_num, n_cores, max_minutes, field_name):
    mpi_per_node = 4*64
    n_nodes = 1 + n_cores // 68
    c_knl = 272 // mpi_per_node
    rho = "81.66"

    input_fname = "plt{0}".format(plt_num)
    dgm_fname = "none"
    intl_fname = "none"
    input_full_fname = "{0}".format(input_fname)
    mt_fname = "$SCRATCH/output/my_answers_knl/amr_si_plt_num_{0}_blocks_{1}_field_{2}.mt".format(plt_num, n_cores, field_name)
    # mt_fname = "none"
    log_fname = "log_cc_plt_num_{0}_cores_{1}_field_{2}_amr.txt".format(plt_num, n_cores, field_name)

    script_template = f"""#!/bin/bash -l

#SBATCH -C knl
#SBATCH -A m3509
#SBATCH -q regular
#SBATCH -N {n_nodes}
#SBATCH -t 00:{max_minutes}:00
#SBATCH -L SCRATCH
#SBATCH -o {log_fname}
#SBATCH -e {log_fname}

srun -n {n_cores} -c {c_knl}  --cpu-bind=cores ~/code/reeber/build_knl_double/examples/amr-connected-components/amr_connected_components_double -f {field_name} --plotfile -i {rho} -n -w {input_full_fname} {mt_fname} {dgm_fname} {intl_fname}\n"""

    with open(file_name, 'w') as sl_file:
        sl_file.write(script_template)


if __name__ == "__main__":
    sbatch_file_name = "queue_test_1024.sh"
    max_minutes = 1
    with open(sbatch_file_name, 'w') as sbatch_file:
        sbatch_file.write('#!/bin/sh\n\n')
        for field_name in ["particle_mass_density", "density"]:
            for plt_num in ["03387"]:
                for n_cores in [1024]:
                    file_name = "run_test_{0}_{1}_{2}_amr.sl".format(plt_num, n_cores, field_name)
                    write_sl(file_name, plt_num, n_cores, max_minutes, field_name)
                    sbatch_file.write("sbatch {0}\n".format(file_name))

