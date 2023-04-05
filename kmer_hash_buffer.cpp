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

    // Create the distributed objects here for buffer
    // Size is n*buffer_size

    std::cout << "a" << std::endl;
    int buffer_size = 10;
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> send_buffer_g(upcxx::new_array<kmer_pair>(buffer_size * upcxx::rank_n()));
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> recv_buffer_g(upcxx::new_array<kmer_pair>(buffer_size * upcxx::rank_n()));

    std::cout << "b" << std::endl;
    // Instantiate the hash table
    HashMap hashmap(hash_table_size, proc_hash_table_size, buffer_size, send_buffer_g, recv_buffer_g);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %d for %d kmers.\n", hash_table_size,
                     n_kmers);
    }

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

    if (run_type == "verbose") {
        BUtil::print("Finished reading kmers.\n");
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "c" << std::endl;
    std::vector<kmer_pair> start_nodes;

    int num_kmers_read = 0;

    while (num_kmers_read < n_kmers) {

        std::cout << num_kmers_read << std::endl;
        kmer_pair kmer = kmers[num_kmers_read];
        
        // Fill the buffers
        bool buffers_full = false;
        while (!buffers_full) {
             buffers_full = hashmap.fill_buffer(kmer);
             num_kmers_read++;  
             if (kmer.backwardExt() == 'F') {
                    start_nodes.push_back(kmer);
                }
        }
        std::cout << num_kmers_read << std::endl;
        
        upcxx::barrier();

        std::cout << "aaa" << std::endl;
    
        // Send messages (one-sided)
        for (int proc_num = 0; proc_num < upcxx::rank_n(); proc_num++) {
            std::cout << "rank " << upcxx::rank_me() << " sending to " << proc_num << std::endl;
            // Get the index of the processor we want to send to
            upcxx::global_ptr<kmer_pair> target_proc_buffer_pointer = recv_buffer_g.fetch(proc_num).wait();
            std::cout << "ccc" << std::endl;
            // Need to define what send_buffer is now
            kmer_pair *send_buffer = hashmap.get_send_buffer(proc_num);
            std::cout << "ddd" << std::endl;
            std::cout << target_proc_buffer_pointer << std::endl;
            for (int j = 0; j < 10; j ++) {
                    std::cout << send_buffer[j].kmer.get() << std::endl;

            }
            //upcxx::rput(&send_buffer, target_proc_buffer_pointer).wait();
            std::cout << "eee" << std::endl;


        }
        std::cout << "bbb" << std::endl;
        upcxx::barrier();
                /*
        // Do local hash inserts

            
        
        */  

    }

/*
    for (auto& kmer : kmers) {

        bool success = hashmap.insert(kmer);
        if (!success) {
            throw std::runtime_error("Error: HashMap is full!");
        }

        
    }
    */

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
