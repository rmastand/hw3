#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {
    size_t my_size;
    size_t size() const noexcept;

    HashMap(size_t size);

    /* 
    "A common way to implement a distributed array in UPC++ is to create a C++ vector of upcxx::global_ptrs 
    that point to arrays allocated in the shared segment on each rank."
    */

    // Each local hash table in the shared memory should be smaller
    size_t shared_table_size = (size / upcxx::rank_n()) + 1

    // Create the distributed objects
    upcxx::dist_object<upcxx::global_ptr<int>> data_g(upcxx::new_array<kmer_pair>(shared_table_size));
    upcxx::dist_object<upcxx::global_ptr<int>> used_g(upcxx::new_array<int>(shared_table_size));

    // Downcast the global pointers
    int *data = data_g->local();
    int *used = used_g->local();

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Helper functions

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t size) {
    my_size = size;
    data.resize(size);
    used.resize(size, 0);
}

bool HashMap::insert(const kmer_pair& kmer) {
    /* 
        Radha says: very confused about what specifically needs to be communicated here.
       I *think* the idea is that if there's no space in the local processes' hash, it needs to either
       send the kmer to the next process, or get ahold of the hash of the next process.
       */
    // Get the hash code
    uint64_t hash = kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
        uint64_t slot = (hash + probe++) % size();
        success = request_slot(slot);
        if (success) {
            write_slot(slot, kmer);
        }
    } while (!success && probe < size());
    return success;
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
        uint64_t slot = (hash + probe++) % size();
        if (slot_used(slot)) {
            val_kmer = read_slot(slot);
            if (val_kmer.kmer == key_kmer) {
                success = true;
            }
        }
    } while (!success && probe < size());
    return success;
}

bool HashMap::slot_used(uint64_t slot) { return used[slot] != 0; }

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data[slot] = kmer; }

kmer_pair HashMap::read_slot(uint64_t slot) { return data[slot]; }

bool HashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0) {
        return false;
    } else {
        used[slot] = 1;
        return true;
    }
}

size_t HashMap::size() const noexcept { return my_size; }
