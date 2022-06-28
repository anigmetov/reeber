#!/usr/bin/env python3
import re

def read_log(fname):
    read_time = local_time  = exch_time = write_time = n_rounds = 0
    try:
        with open(fname, 'r') as f:
            for s in f:
                if s.find("MASTER round") >= 0:
                    n_rounds += 1
                search_result = re.search("read: ([0-9]*) local: ([0-9]*) exchange: ([0-9]*) output: ([0-9]*)", s)
                if search_result:
                    read_time  = int(search_result.group(1))
                    local_time  = int(search_result.group(2))
                    exch_time  = int(search_result.group(3))
                    write_time  = int(search_result.group(4))
    except FileNotFoundError:
        pass
    return [ (local_time + exch_time) / 1000, read_time, local_time, exch_time, write_time, n_rounds // 2 ]


if __name__ == "__main__":
    full_table = {}

    # print(read_log("log_simple_cori_knl_big_plt_num_00632_cores_131072_field_particle_mass_density_amr.txt"))

    plt_nums = [ "01272" ]

    for plt_num in plt_nums:
        for ht in [True, False]:
            for n_cores in [16 * 2048, 32 * 2048, 64 * 2048, 128 * 2048]:
                table_key = (n_cores, plt_num, ht)
                if not ht:
                    log_fname = "log_cc_cori_knl_noht_big_plt_num_{0}_cores_{1}_amr.txt".format(plt_num, n_cores)
                else:
                    log_fname = "log_cc_cori_knl_big_plt_num_{0}_cores_{1}_amr.txt".format(plt_num, n_cores)
                res = read_log(log_fname)
                full_table[table_key] = res


    for ht in [False]:
        for plt_num in plt_nums:
            print("\nHyperthreading: {}, file: {}\n".format(ht, plt_num))
            for n_cores in [16 * 2048, 32 * 2048, 64 * 2048, 128 * 2048]:
                s = ""
                table_key = (n_cores, plt_num, ht)
                for t in full_table[table_key]:
                    s += str(t) + "; "
                print(s)
            print("\n")

