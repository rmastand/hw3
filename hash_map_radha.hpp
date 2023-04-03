#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {

    size_t full_size;
    size_t size() const noexcept;

    size_t local_size;
    size_t local_hash_size() const noexcept;

    // Create the distributed objects
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> *data_g1;
    upcxx::dist_object<upcxx::global_ptr<int>> *used_g1;

     // Create pointers for the data and used objects of the processor that contains the relevant part of the hash map
    upcxx::global_ptr<double> data_pointer;
    upcxx::global_ptr<double> used_pointer;

    // Create objects for the data and used objects of the processor that contains the relevant part of the hash map
    std::vector<kmer_pair> data;
    std::vector<int> used;

    HashMap(size_t full_size, size_t local_size, upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &data_g, upcxx::dist_object<upcxx::global_ptr<int>> &used_g);

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

HashMap::HashMap(size_t full_size, size_t local_size, upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &data_g, upcxx::dist_object<upcxx::global_ptr<int>> &used_g) {
    
    full_size = full_size;
    local_size = local_size;
    data_g1 = &data_g;
    used_g1 = &used_g;

}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    uint64_t probe = 0;
    bool success = false;
    uint64_t slot = (hash + probe) % size();

    // Get the local pointer of the processor that has this slot
    int target_proc_index = slot / local_hash_size();

    // Fetch the pointers for for data and used the target processor
    // What if this proc is full??
    // This should be atomic fetch?
    //data_pointer = data_g1.fetch(target_proc_index).wait();
    used_pointer = used_g1.fetch(target_proc_index).wait();

    // Get the values of data and used in the pointers
    data = upcxx::rget(data_pointer).wait(); // change to atomic!! should this be a future??
    used = upcxx::rget(used_pointer).wait(); // change to atomic!!

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

size_t HashMap::size() const noexcept { return full_size; }

size_t HashMap::local_hash_size() const noexcept { return local_size; }
