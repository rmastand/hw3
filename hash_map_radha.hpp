#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {
    // Size of the full hash table
    size_t my_size_full; 
    size_t full_size() const noexcept; 
    // Size of the shared component of the hash table
    size_t my_size_loc;
    size_t shared_size() const noexcept; 
   
    HashMap(size_t size);

        // Create pointers for the data and used objects of the processor that contains the relevant part of the hash map
    upcxx::global_ptr<double> data_pointer;
    upcxx::global_ptr<double> used_pointer;

    // Create objects for the data and used objects of the processor that contains the relevant part of the hash map
    std::vector<kmer_pair> data;
    std::vector<int> used;


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
    my_size_full = size;

    // Each local hash table in the shared memory should be smaller
    my_size_loc = my_size_full / upcxx::rank_n() + 1;

    // Create the distributed objects
    upcxx::dist_object<upcxx::global_ptr<int>> data_g(upcxx::new_array<kmer_pair>(my_size_loc));
    upcxx::dist_object<upcxx::global_ptr<int>> used_g(upcxx::new_array<int>(my_size_loc));




}

bool HashMap::insert(const kmer_pair& kmer) {
    // Get the hash code
    uint64_t hash = kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    uint64_t slot = (hash + probe) % full_size();

    // Get the local pointer of the processor that has this slot
    int target_proc_index = slot / shared_size();

    // Fetch the pointers for for data and used the target processor
    // What if this proc is full??
    // This should be atomic fetch?
    data_pointer = data_g.fetch(target_proc_index).wait();
    used_pointer = used_g.fetch(target_proc_index).wait();

    // Get the values of data and used in the pointers
    data = upcxx::rget(data_pointer).wait(); // change to atomic!! should this be a future??
    used = upcxx::rget(used_pointer).wait(); // change to atomic!!

    do {
        slot = (hash + probe++) % full_size();  // Increment the slot
        success = request_slot(slot);
        if (success) {
            write_slot(slot, kmer);
        }
    } while (!success && probe < full_size());
    return success;
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    do {
        uint64_t slot = (hash + probe++) % full_size();
        if (slot_used(slot)) {
            val_kmer = read_slot(slot);
            if (val_kmer.kmer == key_kmer) {
                success = true;
            }
        }
    } while (!success && probe < full_size());
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

size_t HashMap::full_size() const noexcept { return my_size_full; }
size_t HashMap::shared_size() const noexcept { return my_size_loc; }
