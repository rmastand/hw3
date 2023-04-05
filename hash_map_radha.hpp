#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>
#include <iostream>

struct HashMap {

    // Create the atomic domain here
    upcxx::atomic_domain<int> ad = upcxx::atomic_domain<int>({upcxx::atomic_op::fetch_add,upcxx::atomic_op::load});

    size_t full_table_size;
    size_t size() const noexcept;

    size_t local_table_size;
    size_t local_size() const noexcept;

    // Create the distributed objects
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> *data_g;
    upcxx::dist_object<upcxx::global_ptr<int>> *used_g;

    HashMap(size_t full_table_size1, size_t local_table_size1, upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &data_g1, upcxx::dist_object<upcxx::global_ptr<int>> &used_g1);
    ~HashMap();    
    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Helper functions

    // Write and read to a logical data slot in the table.
    /*
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
    */
};

HashMap::HashMap(size_t full_table_size1, size_t local_table_size1, upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &data_g1, upcxx::dist_object<upcxx::global_ptr<int>> &used_g1) {
    
    full_table_size = full_table_size1;
    local_table_size = local_table_size1;
    data_g = &data_g1;
    used_g = &used_g1;

}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    bool success = false;
    uint64_t global_slot = hash % size();

   // Useful atomic domain
   //upcxx::atomic_domain<int> ad1({upcxx::atomic_op::fetch_add});

    // Get the index of the processor that has the slot for the hash
    int target_proc_index = global_slot / local_size();

    // We may need to consider multiple target processors as things get filled
    for (int next_proc = 0; next_proc <= upcxx::rank_n(); next_proc++) {

        target_proc_index = (target_proc_index + next_proc) % upcxx::rank_n();

        // Get the pointer for used for that target processor
        upcxx::global_ptr<int> target_proc_used_pointer = used_g->fetch(target_proc_index).wait();

        // Determine the start of where we start looking for empty slots
        uint64_t local_slot;
        if (next_proc == 0) {local_slot =  global_slot % local_size();}
        else {local_slot = 0;}

        for (int probe = 0; (probe + local_slot) < local_size(); probe++) {

            // Iterate through used of the target processor
            int is_slot_full = ad.fetch_add(target_proc_used_pointer+local_slot+probe, 1, std::memory_order_relaxed).wait();
            
            if (is_slot_full == 0) {

                    // Store the kmer
                    upcxx::global_ptr<kmer_pair> target_proc_data_pointer = data_g->fetch(target_proc_index).wait();
                    upcxx::rput(kmer, target_proc_data_pointer+local_slot+probe).wait();
                    success = true;
                    return success;  
                    
                }
            }
        }
        return success;
    }



bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    bool success = false;
    uint64_t global_slot = hash % size();

    // Useful atomic domain
    //upcxx::atomic_domain<int> ad2({upcxx::atomic_op::load});

    // Get the index of the processor that has the slot for the hash
    int target_proc_index = global_slot / local_size();

    // We may need to consider multiple target processors as things get filled
    for (int next_proc = 0; next_proc < upcxx::rank_n() + 1; next_proc++) {

        target_proc_index = (target_proc_index + next_proc) % upcxx::rank_n();

        // Get the pointers for data and used for that target processor
        upcxx::global_ptr<int> target_proc_used_pointer = used_g->fetch(target_proc_index).wait();

        // Determine the start of where we start looking for the hash
        uint64_t local_slot;
        if (next_proc == 0) {local_slot = global_slot % local_size();}
        else {local_slot = 0;}


        for (int probe = 0; (probe + local_slot) < local_size(); probe++) {

            // Iterate through used of the target processor
            int is_slot_full = ad.load(target_proc_used_pointer+local_slot+probe, std::memory_order_relaxed).wait();

            if (is_slot_full > 0) {
                  
                 upcxx::global_ptr<kmer_pair> target_proc_data_pointer = data_g->fetch(target_proc_index).wait();
                val_kmer = upcxx::rget(target_proc_data_pointer+local_slot+probe).wait();
 
                if (val_kmer.kmer == key_kmer) {
                    success = true;
                    return success;      
                } 
                    
                }
            }
        }

        return success;
    }




/*
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
*/

size_t HashMap::size() const noexcept { return full_table_size; }

size_t HashMap::local_size() const noexcept { return local_table_size; }

HashMap::~HashMap() {
    // Need to destroy the atomic domain when the hashmap is destroyed
    ad.destroy();
}
