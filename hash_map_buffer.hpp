#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>
#include <iostream>

#define BUFFER_SIZE 40

struct HashMap {

    size_t full_table_size;
    size_t size() const noexcept;

    size_t local_table_size;
    size_t local_size() const noexcept;

    int* how_many_in_buffer;

    // Downcasted objects
    std::vector<kmer_pair> * send_buff;
    std::vector<kmer_pair> * recv_buff;
    kmer_pair * data_loc;
    int * used_loc;

    // Distributed objects
    upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> *recv_buffer_g;
    upcxx::dist_object<upcxx::global_ptr<kmer_pair>> *data_g;
    upcxx::dist_object<upcxx::global_ptr<int>> *used_g;

    HashMap(size_t full_table_size1, size_t local_table_size1, 
            std::vector<kmer_pair> * send_buff1,
            upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> &recv_buffer_g1,
            upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &data_g1,
            upcxx::dist_object<upcxx::global_ptr<int>> &used_g1);

    ~HashMap();  

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    bool fill_buffer(const kmer_pair& kmer);

    // Functions to deal with the buffers
    bool clear_buffer();
    bool send_all_buffers();
    
    upcxx::future<> send_buffer(upcxx::global_ptr<std::vector<kmer_pair>> remote_dst, int sending_rank, std::vector<kmer_pair> buf);
    upcxx::future<kmer_pair> find_rpc(upcxx::global_ptr<int> remote_dst_used, upcxx::global_ptr<kmer_pair> remote_dst_data, 
                            const pkmer_t &kmer_key_to_find, int slot_to_start, int local_proc_size);
    
    // Helper functions

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t full_table_size1, size_t local_table_size1,
        std::vector<kmer_pair> * send_buff1,
        upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> &recv_buffer_g1,
        upcxx::dist_object<upcxx::global_ptr<kmer_pair>> &data_g1,
        upcxx::dist_object<upcxx::global_ptr<int>> &used_g1) {

    // Constants
    full_table_size = full_table_size1;
    local_table_size = local_table_size1;

    // Buffers
    send_buff = send_buff1;
    recv_buffer_g = &recv_buffer_g1;
    recv_buff = (*recv_buffer_g)->local();

    // Hash table
    data_g = &data_g1;
    used_g = &used_g1;
    data_loc = (*data_g)->local();
    used_loc = (*used_g)->local();

    // Keeps track of how full the buffers are
    how_many_in_buffer = new int[upcxx::rank_n()];
    for (int i = 0; i < upcxx::rank_n(); i++) {
        how_many_in_buffer[i] = 0;
    }

    // Initialize the processor's used to be 0
    for(int i = 0; i < local_table_size; i++) {
      used_loc[i] = 0; 
   }

}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    uint64_t global_slot = hash % size();
    bool success;

    // Get the index of the processor that has the slot for the hash
    int target_proc_index = global_slot / local_size();

    // Add the kmer to the buffer
    send_buff[target_proc_index].push_back(kmer);
    how_many_in_buffer[target_proc_index]++;

    int my_rank = upcxx::rank_me();

    // Send out any filled buffers if buffer size reached
    if (how_many_in_buffer[target_proc_index] == BUFFER_SIZE) {

        // If sending to another processor
        if (target_proc_index != my_rank) {

        // Get the pointer to the remove recv buffer
        upcxx::global_ptr<std::vector<kmer_pair>> target_proc_buffer_pointer = recv_buffer_g->fetch(target_proc_index).wait();

        // Send the buffer
        send_buffer(target_proc_buffer_pointer, my_rank, send_buff[target_proc_index]);//,  send_buff[target_proc_index], sending_rank);
        
        }

        // Place in my own buffer
        else {
            for (auto loc_kmer: send_buff[my_rank]) {
                recv_buff[my_rank].push_back(loc_kmer);
            }
        }

        // Empty the sending buffer
        send_buff[target_proc_index].clear();
        how_many_in_buffer[target_proc_index] = 0;
        }

    // Clear buffer
    bool cleared = clear_buffer();
    return cleared;

}

bool HashMap::send_all_buffers() {
    // cleanup function 
    int my_rank = upcxx::rank_me();

    for (int target_proc_index = 0 ; target_proc_index < upcxx::rank_n(); target_proc_index++) {

        // If sending to another processor
        if (target_proc_index != my_rank) {

            // Get the pointer to the remove recv buffer
            upcxx::global_ptr<std::vector<kmer_pair>> target_proc_buffer_pointer = recv_buffer_g->fetch(target_proc_index).wait();
            
            // Send the buffer
            send_buffer(target_proc_buffer_pointer, my_rank, send_buff[target_proc_index]);
         }

        else {
            for (auto loc_kmer: send_buff[my_rank]) {
                recv_buff[my_rank].push_back(loc_kmer);
            }
        }

        // Empty the sending buffer
        send_buff[target_proc_index].clear();
        how_many_in_buffer[target_proc_index] = 0;
     }

     return true;

}


bool HashMap::clear_buffer() {
    
    // Check if any receive buffers are filled and locally hash
    for (int i = 0 ; i < upcxx::rank_n(); i++) {

        unsigned int recv_buff_size = recv_buff[i].size();

        if (recv_buff_size > 0) {

            // run for loop from 0 to vecSize
            for (unsigned int j = 0; j < recv_buff_size; j++) {

                // Remove the front of the vector
                // Note elements are added to the back of the vector
                kmer_pair working_kmer = recv_buff[i].front();
                recv_buff[i].erase(recv_buff[i].begin());

                uint64_t local_hash = working_kmer.hash();
                uint64_t local_slot = (local_hash % size() ) % local_size();
                uint64_t probe = 0;
                bool success = false;
                do {

                    uint64_t slot = (local_slot + probe++) % local_size();
                    success = request_slot(slot);
                    if (success) {
                        write_slot(slot, working_kmer);
                    }
                } while (!success && probe < local_size());

                // After every local insert, check that the insert was successful
                if (!success) {
                    return false;
                }
            }
    }   
    }
return true;
}


upcxx::future<> HashMap::send_buffer(upcxx::global_ptr<std::vector<kmer_pair>> remote_dst, int sending_rank, std::vector<kmer_pair> buf) {
  
  return upcxx::rpc(remote_dst.where(),
    [](const upcxx::global_ptr<std::vector<kmer_pair>> &dst, const int & sender_rank, std::vector<kmer_pair> buf_in_rpc) {
      
        std::vector<kmer_pair> * dst_recv_buff = dst.local(); 
    
        for (auto loc_kmer: buf_in_rpc) {
            dst_recv_buff[sender_rank].push_back(loc_kmer);
        }

    },
    remote_dst, sending_rank, buf);
}




upcxx::future<kmer_pair> HashMap::find_rpc(upcxx::global_ptr<int> remote_dst_used, upcxx::global_ptr<kmer_pair> remote_dst_data, 
                            const pkmer_t &kmer_key_to_find, int slot_to_start, int local_proc_size) {
  
  return upcxx::rpc(remote_dst_used.where(),
    [](const upcxx::global_ptr<int> dst_used,  const upcxx::global_ptr<kmer_pair> dst_data, 
                           const pkmer_t &kmer_key, int starting_slot, int proc_size) {

        int * dst_used_loc = dst_used.local();
        kmer_pair * dst_data_loc = dst_data.local();
        uint64_t probe = 0;
        bool success = false;
        kmer_pair kmer_val;

        do {
            uint64_t slot = (starting_slot + probe++) % proc_size;
            // See if the slot is used
            int is_slot_used = dst_used_loc[slot];
            if (is_slot_used == 1) {
                // Get the kmer
                kmer_val = dst_data_loc[slot];
                //std::cout << "probing ... " << kmer_val.kmer.get() << std::endl;
                if (kmer_val.kmer == kmer_key) {
                    //std::cout << kmer_val.kmer.get() << std::endl;

                    success = true;

                } 
            }
        } while (!success && probe < proc_size);

         return kmer_val;

    },
    remote_dst_used, remote_dst_data, kmer_key_to_find, slot_to_start, local_proc_size);
}


bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    uint64_t global_slot = hash % size();
    bool success = false;

    // Get the index of the processor that has the slot for the hash
    int target_proc_index = global_slot / local_size();

    int my_rank = upcxx::rank_me();
    uint64_t local_slot = global_slot % local_size();
    uint64_t probe = 0;

    if (target_proc_index != my_rank) {

        // Get the pointers to that target's data and used
        upcxx::global_ptr<int> target_proc_used_pointer = used_g->fetch(target_proc_index).wait();
        upcxx::global_ptr<kmer_pair> target_proc_data_pointer = data_g->fetch(target_proc_index).wait();

        val_kmer = find_rpc(target_proc_used_pointer, target_proc_data_pointer, key_kmer, local_slot, local_size()).wait();

        return true;
        }

    else {
        do {
            uint64_t slot = (local_slot + probe++) % local_size();

            // See if the slot is used
            int is_slot_used = used_loc[slot];
            if (is_slot_used == 1) {

                // Get the kmer
                val_kmer = data_loc[slot];
                if (val_kmer.kmer == key_kmer) {
                    success = true;
                    return success;      
                } 
            }
        } while (!success && probe < local_size());
        return success;

        }
}



bool HashMap::slot_used(uint64_t slot) { return used_loc[slot] != 0; }

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data_loc[slot] = kmer; }

kmer_pair HashMap::read_slot(uint64_t slot) { return data_loc[slot]; }

bool HashMap::request_slot(uint64_t slot) {
    if (used_loc[slot] != 0) {
        return false;
    } else {
        used_loc[slot] = 1;
        return true;
    }
}

size_t HashMap::size() const noexcept { return full_table_size; }

size_t HashMap::local_size() const noexcept { return local_table_size; }

HashMap::~HashMap() {
    delete [] how_many_in_buffer;
}
