// uncomment to switch off sparsification
//#define REEBER_NO_SPARSIFICATION

//#define REEBER_DO_DETAILED_TIMING
#define REEBER_EXTRA_INTEGRAL

#include "reeber-real.h"

// to print nice backtrace on segfault signal
#include <signal.h>
#include <execinfo.h>
#include <cxxabi.h>

#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <diy/io/shared.hpp>
#include <diy/decomposition.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>
#include <error.h>

#include "../../amr-merge-tree/include/fab-block.h"
#include "fab-cc-block.h"
#include "reader-interfaces.h"

#include <output-persistence.h>

#include "../../amr-merge-tree/include/read-npy.h"
#include "../../amr-merge-tree/include/read-hdf5.h"

#include "amr-plot-reader.h"
#include "amr-connected-components-complex.h"

#include <lowfive/vol-base.hpp>
#include <lowfive/vol-metadata.hpp>
#include <lowfive/vol-dist-metadata.hpp>


// block-independent types
using AMRLink = diy::AMRLink;

using BoolVector = diy::RegularDecomposer<diy::DiscreteBounds>::BoolVector;

using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;

#define DIM 3

using FabBlockR = FabBlock<Real, DIM>;

using Block = FabComponentBlock<Real, DIM>;
using Vertex = Block::Vertex;
using Component = Block::Component;
using MaskedBox = Block::MaskedBox;
using GidVector = Block::GidVector;
using GidSet = Block::GidSet;

using TripletMergeTree = Block::TripletMergeTree;
using Neighbor = TripletMergeTree::Neighbor;

herr_t
fail_on_hdf5_error(hid_t stack_id, void*)
{
    H5Eprint(stack_id, stderr);
    fprintf(stderr, "An HDF5 error was detected. Terminating.\n");
    exit(1);
}

struct IsAmrVertexLocal
{
    bool operator()(const Block& b, const Neighbor& from) const
    {
        return from->vertex.gid == b.gid;
    }
};

template<class Real, class LocalFunctor>
struct ComponentDiagramsFunctor
{

    ComponentDiagramsFunctor(Block* b, const LocalFunctor& lf)
            :
            block_(b),
            negate_(b->get_merge_tree().negate()),
            ignore_zero_persistence_(true),
            test_local(lf)
    { }

    void operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        if (!test_local(*block_, from))
        {
            return;
        }

        AmrVertexId current_vertex = from->vertex;

        Real birth_time = from->value;
        Real death_time = through->value;

        if (ignore_zero_persistence_ and birth_time == death_time)
        {
            return;
        }

        AmrVertexId root = block_->vertex_to_deepest_[current_vertex];
        block_->local_diagrams_[root].emplace_back(birth_time, death_time);
    }

    Block* block_;
    const bool negate_;
    const bool ignore_zero_persistence_;
    LocalFunctor test_local;
};

using OutputPairsR = OutputPairs<Block, IsAmrVertexLocal>;

inline bool file_exists(const std::string& s)
{
    std::ifstream ifs(s);
    return ifs.good();

}

inline bool ends_with(const std::string& s, const std::string& suffix)
{
    if (suffix.size() > s.size())
    {
        return false;
    }
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

void read_from_file(std::string infn,
        std::vector<std::string> all_var_names, // HDF5 only: all fields that will be read from plotfile
        int n_mt_vars,                          // sum of first n_mt_vars in all_var_names will be stored in fab of FabBlock,
                                                // for each variable listed in all_var_names FabBlock will have an extra GridRef
        diy::mpi::communicator& world,
        diy::Master& master_reader,
        diy::ContiguousAssigner& assigner,
        diy::MemoryBuffer& header,
        diy::DiscreteBounds& domain,
        bool split,
        int nblocks,
        BoolVector wrap)
{
        read_from_hdf5_file(infn, all_var_names, n_mt_vars, world, nblocks, master_reader, assigner, header, domain);
}

void create_fab_cc_blocks(const diy::mpi::communicator& world, int in_memory, int threads, Real rho, bool absolute,
        bool negate, bool wrap, const diy::FileStorage& storage, diy::Master& master_reader, diy::Master& master,
        const Real cell_volume, const diy::DiscreteBounds& domain)
{
    // copy FabBlocks to FabComponentBlocks
    // in FabTmtConstructor mask will be set and local trees will be computed
    // FabBlock can be safely discarded afterwards

    master_reader.foreach(
            [&master, domain, rho, negate, wrap, absolute, cell_volume](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<AMRLink*>(cp.link());
                AMRLink* new_link = new AMRLink(*l);

                // prepare neighbor box info to save in MaskedBox
                // TODO: refinment vector
                int local_ref = l->refinement()[0];
                int local_lev = l->level();

                master.add(cp.gid(),
                        new Block(b->fab, b->extra_names_, b->extra_fabs_, local_ref, local_lev, domain,
                                l->bounds(),
                                l->core(), cp.gid(),
                                new_link, rho, negate, wrap, absolute, cell_volume),
                        new_link);

            });
}

void write_tree_blocks(const diy::mpi::communicator& world, bool split, const std::string& output_filename,
        diy::Master& master)
{
    if (output_filename != "none")
    {
        if (!split)
        {
            diy::io::write_blocks(output_filename, world, master);
        } else
        {
            diy::io::split::write_blocks(output_filename, world, master);
        }
    }
}

extern "C" void consumer_f(int argc, char** argv, MPI_Comm& comm, MPI_Comm& intercomm)
{
    LowFive::create_logger("trace");
    // setup low five
    LowFive::LocationPattern only_one {"plt00002.h5", "*"};
    H5Eset_auto(H5E_DEFAULT, fail_on_hdf5_error, NULL);
    LowFive::DistMetadataVOL vol(comm, intercomm);
    vol.set_intercomm("plt00002.h5", "*", 0);
    vol.set_memory("*", "*");

    diy::mpi::communicator world(comm);

    int nblocks = world.size();
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    int threads = 1;
    std::string profile_path;
    std::string log_level = "info";

    dlog::add_stream(std::cerr, dlog::severity(log_level))
            << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level()
            << dlog::color_post() >> dlog::flush();

    // threshold
    Real rho = 81.66;
    Real absolute_rho;
    int min_cells = 10;

    std::string integral_fields = "/level_0/density";
    std::string function_fields = "/level_0/density";

    int n_runs = 1;

    using namespace opts;

    LOG_SEV_IF(world.rank() == 0, info) << "Started, n_runs = " << n_runs << ", size of Real = " << sizeof(Real);

    opts::Options ops(argc, argv);
    ops
            >> Option('b', "blocks", nblocks, "number of blocks to use")
            >> Option('m', "memory", in_memory, "maximum blocks to store in memory")
            >> Option('j', "jobs", threads, "threads to use during the computation")
            >> Option('s', "storage", prefix, "storage prefix")
            >> Option('i', "rho", rho, "iso threshold")
            >> Option('x', "mincells", min_cells, "minimal number of cells to output halo")
            >> Option('f', "function_fields", function_fields, "fields to add for merge tree, separated with , ")
            >> Option(      "integral_fields", integral_fields, "fields to integrate separated with , ")
            >> Option('r', "runs", n_runs, "number of runs")
            >> Option('p', "profile", profile_path, "path to keep the execution profile")
            >> Option('l', "log", log_level, "log level");

    bool absolute =
            ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean)");
    bool negate = ops >> opts::Present('n', "negate", "sweep superlevel sets");
    // ignored for now, wrap is always assumed
    bool wrap = ops >> opts::Present('w', "wrap", "wrap");
    bool split = ops >> opts::Present("split", "use split IO");

    BoolVector wrap_vec { wrap, wrap, wrap };

    bool print_stats = ops >> opts::Present("stats", "print statistics");
    std::string input_filename = "plt00002.h5";
    std::string output_filename = "out_plt00002.h5";
    std::string output_diagrams_filename, output_integral_filename;

    std::vector<std::string> all_var_names;
    std::vector<std::string> integral_var_names;
    int n_mt_vars;
    Real cell_volume = 1.0;

    // we read hdf5 only
    std::vector<std::string> function_var_names = split_by_delim(function_fields, ',');  //{"particle_mass_density", "density"};
    n_mt_vars = function_var_names.size();
    if (integral_fields.empty())
    {
        all_var_names = function_var_names;
    } else
    {
        all_var_names = function_var_names;
        integral_var_names = split_by_delim(integral_fields, ',');  //{"particle_mass_density", "density"};
        for(std::string i_name : integral_var_names)
        {
            if (std::find(function_var_names.begin(), function_var_names.end(), i_name) == function_var_names.end())
            {
                all_var_names.push_back(i_name);
            }
        }
    }

    LOG_SEV_IF(world.rank() == 0, info) << "Reading fields: " << container_to_string(all_var_names) << ", fields to sum = " << n_mt_vars;
    dlog::flush();

    bool write_diag = false;
    bool write_integral = false;

    diy::FileStorage storage(prefix);

    diy::Master master_reader(world, 1, in_memory, &FabBlockR::create, &FabBlockR::destroy);
    diy::ContiguousAssigner assigner(world.size(), nblocks);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain(DIM);

    dlog::Timer timer;
    dlog::Timer timer_all;
    LOG_SEV_IF(world.rank() == 0, info) << "Starting computation, input_filename = " << input_filename << ", nblocks = "
                                        << nblocks
                                        << ", rho = " << rho;
    dlog::flush();

   read_from_file(input_filename, all_var_names, n_mt_vars, world, master_reader, assigner, header, domain, split, nblocks, wrap_vec);

    auto time_to_read_data = timer.elapsed();
    dlog::flush();

    LOG_SEV_IF(world.rank() == 0, info) << "Data read, local size = " << master_reader.size();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:       " << dlog::clock_to_string(timer.elapsed());
    dlog::flush();

    timer.restart();

    for(int n_run = 0; n_run < n_runs; ++n_run)
    {
        world.barrier();
        timer.restart();
        timer_all.restart();

        diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save,
            &Block::load);

        create_fab_cc_blocks(world, in_memory, threads, rho, absolute, negate, wrap, storage, master_reader, master, cell_volume, domain);

        auto time_for_local_computation = timer.elapsed();



        Real mean = std::numeric_limits<Real>::min();
        Real total_sum;

        if (absolute)
        {
            absolute_rho = rho;
            LOG_SEV_IF(world.rank() == 0, info) << "Time to compute local trees and components:  "
                                                << dlog::clock_to_string(timer.elapsed());
        } else
        {
            LOG_SEV_IF(world.rank() == 0, info) << "Time to construct FabComponentBlocks: "
                                                << dlog::clock_to_string(timer.elapsed());
            dlog::flush();
            timer.restart();

            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
                cp.collectives()->clear();
                cp.all_reduce(b->sum_ * b->scaling_factor(), std::plus<Real>());
                cp.all_reduce(static_cast<Real>(b->n_unmasked_) * b->scaling_factor(), std::plus<Real>());
            });

            master.exchange();

            const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

            total_sum = proxy.get<Real>();
            Real total_unmasked = proxy.get<Real>();

            mean = total_sum / total_unmasked;

            absolute_rho = rho * mean;                                            // now rho contains absolute threshold

            LOG_SEV_IF(world.rank() == 0, info) << "Total sum = " << total_sum << ", total_unmasked = "
                                                << total_unmasked;

            LOG_SEV_IF(world.rank() == 0, info) << "Average = " << mean << ", rho = " << rho
                                                << ", absolute_rho = " << absolute_rho
                                                << ", time to compute average: "
                                                << dlog::clock_to_string(timer.elapsed());

            if (mean < 0 or std::isnan(mean) or std::isinf(mean) or mean > 1e+40)
            {
                LOG_SEV_IF(world.rank() == 0, error) << "Bad average = " << mean << ", do not proceed";
                return;
            }

            time_for_local_computation += timer.elapsed();
            dlog::flush();
            timer.restart();

            long int local_active = 0;
            master.foreach([absolute_rho, &local_active](Block* b, const diy::Master::ProxyWithLink& cp)
            {
                AMRLink* l = static_cast<AMRLink*>(cp.link());
                b->init(absolute_rho, l, true);
                cp.collectives()->clear();
                local_active += b->n_active_;
            });

            dlog::flush();


            LOG_SEV_IF(world.rank() == 0, info) << "Time to initialize FabComponentBlocks (low vertices, local trees, components, outgoing edges): " << timer.elapsed();
            time_for_local_computation += timer.elapsed();
        }

        dlog::flush();

        timer.restart();

        int global_n_undone = 1;

        master.foreach(&send_edges_to_neighbors_cc<Real, DIM>);
        master.exchange();
        master.foreach(&delete_low_edges_cc<Real, DIM>);

        LOG_SEV_IF(world.rank() == 0, info) << "edges symmetrized, time elapsed " << timer.elapsed();
        auto time_for_communication = timer.elapsed();
        dlog::flush();
        timer.restart();

        int rounds = 0;
        while(global_n_undone)
        {
            rounds++;
            master.foreach(&amr_cc_send<Real, DIM>);
            master.exchange();
            master.foreach(&amr_cc_receive<Real, DIM>);
            LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", get OK";
            dlog::flush();
            master.exchange();
            global_n_undone = master.proxy(master.loaded_block()).read<int>();

            LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", global_n_undone = "
                                                << global_n_undone;

            if (print_stats)
            {
                int local_n_undone = 0;
                master.foreach(
                        [&local_n_undone](Block* b, const diy::Master::ProxyWithLink& cp) {
                            local_n_undone += (b->done_ != 1);
                        });

                LOG_SEV(info) << "STAT MASTER round " << rounds << ", rank = " << world.rank() << ", local_n_undone = "
                              << local_n_undone;
            }
            dlog::flush();
        }

        // here merge tree computation stops
        // sync to measure runtime
        world.barrier();
        auto time_total_computation = timer_all.elapsed();

        //    fmt::print("world.rank = {}, time for exchange = {}\n", world.rank(), dlog::clock_to_string(timer.elapsed()));

        LOG_SEV_IF(world.rank() == 0, info) << "Time for exchange:  " << dlog::clock_to_string(timer.elapsed());
        LOG_SEV_IF(world.rank() == 0, info) << "Total time for computation:  " << time_total_computation;
        time_for_communication += timer.elapsed();
        dlog::flush();
        timer.restart();

        // save the result
        write_tree_blocks(world, split, output_filename, master);

        LOG_SEV_IF(world.rank() == 0, info) << "Time to write tree:  " << dlog::clock_to_string(timer.elapsed());
        auto time_for_output = timer.elapsed();
        dlog::flush();
        timer.restart();

        bool verbose = false;

        if (write_diag)
        {
            bool ignore_zero_persistence = true;
            OutputPairsR::ExtraInfo extra(output_diagrams_filename, verbose, world);
            IsAmrVertexLocal test_local;
            master.foreach(
                    [&extra, &test_local, ignore_zero_persistence, absolute_rho](Block* b,
                            const diy::Master::ProxyWithLink& cp)
                    {
                        output_persistence(b, cp, &extra, test_local, absolute_rho, ignore_zero_persistence);
                        dlog::flush();
                    });
        }

        LOG_SEV_IF(world.rank() == 0, info) << "Time to write diagrams:  " << dlog::clock_to_string(timer.elapsed());
        time_for_output += timer.elapsed();
        dlog::flush();
        timer.restart();

#ifdef REEBER_EXTRA_INTEGRAL
        if (write_integral)
        {
            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp)
            {
                b->compute_final_connected_components();
                b->compute_local_integral();
                b->multiply_integral_by_cell_volume();
            });

            LOG_SEV_IF(world.rank() == 0, info) << "Local integrals computed\n";
            dlog::flush();
            world.barrier();

            diy::io::SharedOutFile ofs(output_integral_filename, world);

            master.foreach(
                    [&world, &ofs, domain, min_cells, integral_var_names](Block* b, const diy::Master::ProxyWithLink& cp)
                    {
                        diy::Point<int, 3> domain_shape;
                        for(int i = 0; i < 3; ++i)
                        {
                            domain_shape[i] = domain.max[i] - domain.min[i] + 1;
                        }

                        diy::GridRef<void*, 3> domain_box(nullptr, domain_shape, /* c_order = */ false);

                        // local integral already stores number of vertices (set in init)
                        // so we add it here to the list of fields just to print it

                        std::vector<std::string> integral_vars = b->extra_names_;

                        integral_vars.insert(integral_vars.begin(), std::string("total_mass"));
                        integral_vars.insert(integral_vars.begin(), std::string("n_vertices"));
                        integral_vars.insert(integral_vars.begin(), std::string("n_cells"));

                        LOG_SEV_IF(world.rank() == 0, debug) << "integral_vars:  " << container_to_string(integral_vars);

                        bool print_header = false;
                        if (print_header)
                        {
                            std::string integral_header = "# id x y z  ";
                            for(auto s : integral_vars)
                            {
                                integral_header += s;
                                integral_header += " ";
                            }
                            integral_header += "\n";

                            fmt::print(ofs, integral_header);
                        }

                        for(const auto& root_values_pair : b->local_integral_)
                        {
                            AmrVertexId root = root_values_pair.first;
                            if (root.gid != b->gid)
                                continue;

                            auto& values = root_values_pair.second;

                            Real n_cells = values.at("n_cells");

                            if (n_cells < min_cells)
                                continue;

                            auto root_position = coarsen_point(b->local_.global_position(root), b->refinement(), 1);

                            fmt::print(ofs, "{} {} {}\n",
                                    domain_box.index(root_position),
                                    root_position,
                                    b->pretty_integral(root, integral_vars));
                        }
                    });

            LOG_SEV_IF(world.rank() == 0, info) << "Time to compute and write integral:  "
                                                << dlog::clock_to_string(timer.elapsed());
            time_for_output += timer.elapsed();
            dlog::flush();
            timer.restart();
        }
#endif
        dlog::flush();

        world.barrier();

        std::string final_timings = fmt::format("run: {} read: {} local: {} exchange: {} output: {} total: {}\n",
                n_run, time_to_read_data, time_for_local_computation, time_for_communication, time_for_output,
                time_total_computation);
        LOG_SEV_IF(world.rank() == 0, info) << final_timings;

        dlog::flush();

        if (print_stats)
        {
            std::size_t max_n_gids = 0;
            std::set<int> gids;
            master.foreach([&max_n_gids, &gids](Block* b, const diy::Master::ProxyWithLink& cp) {
                std::set<int> block_gids;
                for(const Component& c : b->components_)
                {
                    gids.insert(c.current_neighbors().begin(), c.current_neighbors().end());
                    block_gids.insert(c.current_neighbors().begin(), c.current_neighbors().end());
                }
                max_n_gids = std::max(max_n_gids, static_cast<decltype(max_n_gids)>(block_gids.size()));
            });

            LOG_SEV_IF(max_n_gids > 0, info) << "STAT max_n_gids[" << world.rank() << "] = " << max_n_gids;
            LOG_SEV_IF(max_n_gids > 0, info) << "STAT total_n_gids[" << world.rank() << "] = " << gids.size();
            dlog::flush();
            world.barrier();

            LOG_SEV_IF(world.rank() == 0, info) << "STAT sizes = np.array(sizes)";
            LOG_SEV_IF(world.rank() == 0, info) << "STAT max_n_gids = np.array(max_n_gids)";
            LOG_SEV_IF(world.rank() == 0, info) << "STAT total_n_gids = np.array(total_n_gids)";
            LOG_SEV_IF(world.rank() == 0, info) << "STAT hist_array = sizes";
            LOG_SEV_IF(world.rank() == 0, info) << "STAT plt.hist(hist_array, bins = 'auto')";
            LOG_SEV_IF(world.rank() == 0, info) << "STAT plt.title('{} cores'.format(n_cores))";
            LOG_SEV_IF(world.rank() == 0, info) << "STAT plt.show()";
        }

        size_t local_n_active = 0;
        size_t local_n_components = 0;
        size_t local_n_blocks = 0;

        if (true)
        {
            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
                size_t n_low = b->n_low_;
                size_t n_active = b->n_active_;
                size_t n_masked = b->n_masked_;
                cp.collectives()->clear();
                cp.all_reduce(n_low, std::plus<size_t>());
                cp.all_reduce(n_active, std::plus<size_t>());
                cp.all_reduce(n_masked, std::plus<size_t>());
            });

            master.exchange();

            world.barrier();

            const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

            size_t total_n_low = proxy.get<size_t>();
            size_t total_n_active = proxy.get<size_t>();
            size_t total_n_masked = proxy.get<size_t>();

            master.foreach([&local_n_active, &local_n_components, &local_n_blocks](Block* b,
                    const diy::Master::ProxyWithLink& cp) {
                local_n_active += b->n_active_;
                local_n_components += b->components_.size();
                local_n_blocks += 1;
            });

            LOG_SEV_IF(world.rank() == 0, info) << "Total_n_low = " << total_n_low << ", total_n_active = "
                                                << total_n_active << ", total_n_masked = "
                                                << total_n_masked;
            dlog::flush();
            timer.restart();
        }

    } // loop over runs

    return;
}
