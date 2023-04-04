#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>

#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"

#include "butil.hpp"
#include <iostream>

int main(int argc, char** argv) {
    upcxx::init();

    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = std::string(argv[1]);
    std::string run_type = "";

    if (argc >= 3) {
        run_type = std::string(argv[2]);
    }

    std::string test_prefix = "test";
    if (run_type == "test" && argc >= 4) {
        test_prefix = std::string(argv[3]);
    }

    int ks = kmer_size(kmer_fname);

    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) +
                                 "-mers.  Modify packing.hpp and recompile.");
    }

    size_t n_kmers = line_count(kmer_fname);

    // Load factor of 0.5
    size_t hash_table_size = n_kmers * (1.0 / 0.5);

    // Size of each processor's hash table
    size_t proc_hash_table_size = hash_table_size / upcxx::rank_n() + 1;
    // if (upcxx::rank_me() == 0) std::cout << hash_table_size <<  " " << proc_hash_table_size << std::endl;
 

    // Create the distributed objects here for data and used
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> data_g(upcxx::new_array<kmer_pair>(proc_hash_table_size));
    upcxx::dist_object<upcxx::global_ptr<int>> used_g(upcxx::new_array<int>(proc_hash_table_size));

    // Initialize the processor's used to be 0
    int *used = used_g->local();
    for(int i = 0; i < proc_hash_table_size; i++) {
      used[i] = 0; 
   }

    // Instantiate the hash table
    HashMap hashmap(hash_table_size, proc_hash_table_size, data_g, used_g);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %d for %d kmers.\n", hash_table_size,
                     n_kmers);
    }

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

    if (run_type == "verbose") {
        BUtil::print("Finished reading kmers.\n");
    }

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;


    for (auto& kmer : kmers) {
        //BUtil::print(kmer.kmer.get());
        //BUtil::print("\n");
        //std::cout << "inserting kmer " << kmer.kmer.get() << " with forward extension " << kmer.forwardExt() << " and backwrd extension " << kmer.backwardExt() << std::endl;
        bool success = hashmap.insert(kmer);
        if (!success) {
            throw std::runtime_error("Error: HashMap is full!");
        }

        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
            BUtil::print("Found a start node 'n");
            BUtil::print(kmer.kmer.get());
            BUtil::print("\n");

        }
    }

    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    double insert_time = std::chrono::duration<double>(end_insert - start).count();
    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf\n", insert_time);
    }
    upcxx::barrier();

    auto start_read = std::chrono::high_resolution_clock::now();

    std::list<std::list<kmer_pair>> contigs;

    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair kmer;
            std::cout << "looking for kmer " << contig.back().next_kmer().get() << " with forward extension " << contig.back().forwardExt() << " and backwrd extension " << contig.back().backwardExt() << std::endl;
            std::cout << "and " << kmer.kmer.get() << std::endl;

            bool success = hashmap.find(contig.back().next_kmer(), kmer);
            if (!success) {
                throw std::runtime_error("Error: k-mer not found in hashmap.");
            }
            contig.push_back(kmer);
        }
        contigs.push_back(contig);
    }

 

    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> read = end_read - start_read;
    std::chrono::duration<double> insert = end_insert - start;
    std::chrono::duration<double> total = end - start;

    int numKmers = std::accumulate(
        contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair>& contig) { return sum + contig.size(); });

    if (run_type != "test") {
        BUtil::print("Assembled in %lf total\n", total.count());
    }

    if (run_type == "verbose") {
        printf("Rank %d reconstructed %d contigs with %d nodes from %d start nodes."
               " (%lf read, %lf insert, %lf total)\n",
               upcxx::rank_me(), contigs.size(), numKmers, start_nodes.size(), read.count(),
               insert.count(), total.count());
    }

    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto& contig : contigs) {
            fout << extract_contig(contig) << std::endl;
        }
        fout.close();
    }

    upcxx::finalize();
    return 0;
}
