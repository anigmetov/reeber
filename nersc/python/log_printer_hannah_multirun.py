#!/usr/bin/env python
import re
import numpy as np

hannah_plt_files = [ "/project/projectdirs/m2848/hannah_Nyx/nbody/1species_vog_256_refined2level_8crit/plotfiles/plt00719" ]

hannah_plt_nums = [ "00719" ]

rhos = [ "81.66" ]
n_coreses = [761, 1522, 3045, 6089, 6089*2]
n_coreses = [2048]

def read_log(fname):
    read_time  = local_time  = exch_time = write_time  = n_rounds = 0
    result = []
    try:
        with open(fname, 'r') as f:
            for s in f:
                if s.find("MASTER round") >= 0 and s.find("get OK") > 0:
                    n_rounds += 1
                search_result = re.search("run: ([0-9]*) read: ([0-9]*) local: ([0-9]*) exchange: ([0-9]*) output: ([0-9]*)", s)
                if search_result:
                    n_run  = int(search_result.group(1))
                    read_time  = int(search_result.group(2))
                    local_time  = int(search_result.group(3))
                    exch_time  = int(search_result.group(4))
                    write_time  = int(search_result.group(5))
                    result.append([ (local_time + exch_time) / 1000, read_time, local_time, exch_time, write_time, n_rounds, n_run ])
                    n_rounds = 0
    except:
        pass
    return result

if __name__ == "__main__":
    full_table = {}

# log_cc_cori_knl_hannah_plt_num_01310_cores_2048_run_2_rho_10_amr.txt
    n_runs = [ 50 ]

    for plt_num in hannah_plt_nums:
        for n_cores in n_coreses:
            for n_run in n_runs:
                for rho in rhos:
                    table_key = (n_cores, plt_num, n_run, rho)
                    log_fname = "log_cc_cori_knl_hannah_plt_num_{0}_cores_{1}_run_{2}_rho_{3}_amr.txt".format(plt_num, n_cores, n_run, rho)
                    res = read_log(log_fname)
                    full_table[table_key] = res


    # for plt_num in hannah_plt_nums:
    #     for rho in rhos:
    #         print("\nTiming for Simple, Hannah, rho = {}, plt_num = {}, average over 9 runs\n".format(rho, plt_num))
    #         for n_cores in n_coreses:
    #             local_avg = 0.0
    #             exch_avg = 0.0
    #             local_cnt = 0
    #             for n_run in n_runs:
    #                 try:
    #                     s = "cores: {}, run: {} | ".format(n_cores, n_run)
    #                     table_key = (n_cores, plt_num, n_run, rho)
    #                     table_entries = full_table[table_key]
    #                     first_entry = True
    #                     for table_entry in table_entries:
    #                         if first_entry:
    #                             first_entry = False
    #                             continue
    #                         s = "cores: {}, run: {} | ".format(n_cores, n_run)
    #                         if table_entry[2] != 0:
    #                             local_avg += table_entry[2]
    #                             exch_avg += table_entry[3]
    #                             local_cnt += 1
    #                         for t in table_entry:
    #                             s += str(t) + "; "
    #                         # print(s)
    #                 except KeyError:
    #                     pass 
    #                 # print("HERE local_cnt = {}".format(local_cnt))
    #             if local_cnt > 0:
    #                 local_avg  /= local_cnt
    #                 exch_avg /= local_cnt
    #                 # print("cores: {0}, local_average = {1:.2f}, exchange_average = {2:.2f}, total average = {3:.2f}".format(
    #                 #     n_cores, local_avg / 1000, exch_avg / 1000, (local_avg + exch_avg) / 1000) )
    #                 print("{0: 5d};{1:.2f}".format(n_cores, (local_avg + exch_avg) / 1000))



    for plt_num in hannah_plt_nums:
        for rho in rhos:
            print("\nTiming for Components, Hannah, rho = {}, plt_num = {}, average over 10 runs\n".format(rho, plt_num))
            for n_cores in n_coreses:
                local_avg = 0.0
                exch_avg = 0.0
                local_cnt = 0
                total_max = 0
                total_min = 10000000
                for n_run in n_runs:
                    try:
                        s = "cores: {}, run: {} | ".format(n_cores, n_run)
                        table_key = (n_cores, plt_num, n_run, rho)
                        table_entries = full_table[table_key]
                        values = []
                        for table_entry in table_entries:
                            s = "cores: {}, run: {} | ".format(n_cores, n_run)
                            if table_entry[2] != 0:
                                local_avg += table_entry[2] / 1000.0
                                exch_avg += table_entry[3] / 1000.0
                                total = (table_entry[2] + table_entry[3]) / 1000.0
                                total_max = max(total, total_max)
                                total_min = min(total, total_min)
                                values.append(total)
                                local_cnt += 1
                            for t in table_entry:
                                s += str(t) + "; "
                            # print(s)
                    except KeyError:
                        pass 
                    # print("HERE local_cnt = {}".format(local_cnt))
                if local_cnt > 0:
                    local_avg  /= local_cnt
                    exch_avg /= local_cnt
                    print("Hee", values)
                    error = np.std(values)
                    # print("cores: {0}, local_average = {1:.2f}, exchange_average = {2:.2f}, total average = {3:.2f}".format(
                    #     n_cores, local_avg / 1000, exch_avg / 1000, (local_avg + exch_avg) / 1000) )
                    # print("{0: 5d};{1:.2f}".format(n_cores, (local_avg + exch_avg) / 1000))
                    print("{0: 5d};{1:.2f};{2:.2f}".format(n_cores,
                        (local_avg + exch_avg) , error))

    # for plt_num in hannah_plt_nums:
    #     for rho in rhos:
    #         print("\nTiming for Simple, Hannah, rho = {}, plt_num = {}, max time\n".format(rho, plt_num))
    #         for n_cores in n_coreses:
    #             total_max = 0
    #             for n_run in n_runs:
    #                 try:
    #                     s = "cores: {}, run: {} | ".format(n_cores, n_run)
    #                     table_key = (n_cores, plt_num, n_run, rho)
    #                     table_entries = full_table[table_key]
    #                     for table_entry in table_entries:
    #                         s = "cores: {}, run: {} | ".format(n_cores, n_run)
    #                         if table_entry[2] != 0:
    #                             if total_max == 0 or table_entry[2] + table_entry[3] > total_max:
    #                                 total_max = table_entry[2] + table_entry[3]
    #                         for t in table_entry:
    #                             s += str(t) + "; "
    #                         # print(s)
    #                 except KeyError:
    #                     pass 
    #                 # print("HERE local_cnt = {}".format(local_cnt))
    #             if local_cnt > 0:
    #                 # print("cores: {0}, local_average = {1:.2f}, exchange_average = {2:.2f}, total average = {3:.2f}".format(
    #                 #     n_cores, local_avg / 1000, exch_avg / 1000, (local_avg + exch_avg) / 1000) )
    #                 print(";{0:.2f}".format(total_max / 1000))
                               
