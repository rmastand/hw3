#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>
#include <iostream>

struct HashMap {

    // Create the atomic domain here
    upcxx::atomic_domain<int> ad = upcxx::atomic_domain<int>({upcxx::atomic_op::compare_exchange,upcxx::atomic_op::load});

    size_t full_table_size;
    size_t size() const noexcept;

    size_t local_table_size;
    size_t local_size() const noexcept;

    std::vector<kmer_pair> data;
    std::vector<int> used;

    int buffer_size;
    int* how_many_in_buffer;
    // Downcast buffers
    kmer_pair *send_buff;
    kmer_pair *recv_buff;



    // Create the distributed objects
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> *send_buffer_g;
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> *recv_buffer_g;

    HashMap(size_t full_table_size1, size_t local_table_size1, int buffer_size, upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &send_buffer_g1,  upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &recv_buffer_g1);
    ~HashMap();  

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    bool fill_buffer(const kmer_pair& kmer);
    bool local_inserts();
    kmer_pair * get_send_buffer(int target_proc);

    // Helper functions

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t full_table_size1, size_t local_table_size1, int buffer_size1, upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &send_buffer_g1, upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &recv_buffer_g1) {
    
    full_table_size = full_table_size1;
    local_table_size = local_table_size1;
    buffer_size = buffer_size1;

    send_buffer_g = &send_buffer_g1;
    recv_buffer_g = &recv_buffer_g1;

    // Keeps track of how full the buffers are
    how_many_in_buffer = new int[upcxx::rank_n()];
    for (int i = 0; i < upcxx::rank_n(); i++) {
        how_many_in_buffer[i] = 0;
    }

    send_buff = (*send_buffer_g)->local();
    recv_buff = (*recv_buffer_g)->local();

    data.resize(local_table_size);
    used.resize(local_table_size, 0);

}

bool HashMap::fill_buffer(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    uint64_t global_slot = hash % size();

    // Get the index of the processor that has the slot for the hash
    int target_proc_index = global_slot / local_size();

    // See where to insert the kmer in the buffer
    int insert_index = buffer_size*target_proc_index + how_many_in_buffer[target_proc_index];
    send_buff[insert_index] = kmer;
    how_many_in_buffer[target_proc_index]++;

    if (how_many_in_buffer[target_proc_index] == buffer_size) {
        return true;
    }

    return false;

}

bool HashMap::local_inserts() {

    bool success;

    for (int ii = 0; ii < buffer_size*upcxx::rank_n(); ii++) {

        kmer_pair kmer = recv_buff[ii];

        int probe = 0;
        uint64_t hash = kmer.hash();
        uint64_t global_slot = hash % size();
        bool success = false;

         // Determine the start of where we start looking for empty slots
        uint64_t local_slot =  global_slot % local_size();

        do {
            uint64_t slot = (local_slot + probe++) % local_size();;
            success = request_slot(slot);
            if (success) {
                write_slot(slot, kmer);
            }
                } while (!success && probe < size());
        }

        return success;

    }




/*

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    bool success = false;
    uint64_t global_slot = hash % size();

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
*/     


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

size_t HashMap::size() const noexcept { return full_table_size; }

size_t HashMap::local_size() const noexcept { return local_table_size; }

kmer_pair* HashMap::get_send_buffer(int target_proc) { 

    kmer_pair * tmp[buffer_size];

    for (int i = 0; i < buffer_size; i++) {
        std::cout <<  "k" << std::endl;
        tmp[i] = send_buff[buffer_size*target_proc + i];
        std::cout <<  send_buff[buffer_size*target_proc + i].kmer.get() << std::endl;
    }

    return *tmp; 
}

HashMap::~HashMap() {
    // Need to destroy the atomic domain when the hashmap is destroyed
    ad.destroy();
    delete [] how_many_in_buffer;
}
