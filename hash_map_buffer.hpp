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

    size_t buffer_size;
    int* how_many_in_buffer;

     // Downcasted buffers
    //std::vector<kmer_pair>(upcxx::rank_n()) send_buff;
    std::vector<kmer_pair> * send_buff;
    //std::vector<kmer_pair>(upcxx::rank_n()) recv_buff;

    // Create the distributed objects
    //upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> *send_buffer_g;
    upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> *recv_buffer_g;

    HashMap(size_t full_table_size1, size_t local_table_size1, size_t buffer_size, 
            std::vector<kmer_pair> * send_buff1,
       // upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> &send_buffer_g1, 
        upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> &recv_buffer_g1);
    ~HashMap();  

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    bool fill_buffer(const kmer_pair& kmer);
    bool local_inserts();
    
    upcxx::future<> send_buffer(upcxx::global_ptr<std::vector<kmer_pair>> remote_dst, int sending_rank, std::vector<kmer_pair> buf);
   // upcxx::future<> send_buffer(upcxx::global_ptr<std::vector<kmer_pair>> receiving_proc, std::vector<kmer_pair> buf, int sending_proc);
    // Helper functions

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t full_table_size1, size_t local_table_size1, size_t buffer_size1, 
            std::vector<kmer_pair> * send_buff1,
            //upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> &send_buffer_g1,
             upcxx::dist_object<upcxx::global_ptr<std::vector<kmer_pair>>> &recv_buffer_g1) {
    full_table_size = full_table_size1;
    local_table_size = local_table_size1;
    buffer_size = buffer_size1;

    //send_buff = (*send_buffer_g)->local();
    //send_buff = send_buffer_g->local();
    //std::cout << "3" << std::endl;
    //recv_buff = (*recv_buffer_g)->local();
    send_buff = send_buff1;
    //send_buffer_g = &send_buffer_g1;
    recv_buffer_g = &recv_buffer_g1;
    // Keeps track of how full the buffers are
    how_many_in_buffer = new int[upcxx::rank_n()];
    for (int i = 0; i < upcxx::rank_n(); i++) {
        how_many_in_buffer[i] = 0;
    }
    // local copies of data and used
    data.resize(local_table_size);
    used.resize(local_table_size, 0);

}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    uint64_t global_slot = hash % size();

    // Get the index of the processor that has the slot for the hash
    int target_proc_index = global_slot / local_size();

    // Add the kmer to the buffer
    send_buff[target_proc_index].push_back(kmer);
    how_many_in_buffer[target_proc_index]++;

    // Send out any filled buffers
    if (how_many_in_buffer[target_proc_index] == buffer_size) {
        // Now the buffer is full, so we need to do a rpc

        std::cout << upcxx::rank_me() << " has full buffer and is sending to " << target_proc_index << std::endl;

        // Get the pointer to the remove recv buffer
        upcxx::global_ptr<std::vector<kmer_pair>> target_proc_buffer_pointer = recv_buffer_g->fetch(target_proc_index).wait();
        // Get the sending processor rank so I can tell the recv buffer where to insert
        int sending_rank = upcxx::rank_me();
        // Send the buffer
        send_buffer(target_proc_buffer_pointer, sending_rank, send_buff[target_proc_index]);//,  send_buff[target_proc_index], sending_rank);
        
        // Empty the sending buffer
        send_buff[target_proc_index].clear();
        how_many_in_buffer[target_proc_index] = 0;

    }

    /*
    // Check if any receive buffers are filled and locally hash
    for (int i = 0 ; i < upcxx::rank_me(); i++) {

        unsigned int recv_buff_size = recv_buff[i].size();

        if (recv_buff_size > 0) {

            // run for loop from 0 to vecSize
            for (unsigned int j = 0; j < recv_buff_size; j++) {

                // Remove the front of the vector
                // Note elements are added to the back of the vector
                kmer_pair working_kmer = recv_buff[i].front();
                recv_buff[i].erase(recv_buff[i].begin());

                uint64_t local_hash = working_kmer.hash();
                uint64_t local_slot = (hash % size() ) % local_size();
                uint64_t probe = 0;
                bool success = false;
                do {
                    uint64_t slot = (local_slot + probe++) % local_size();
                    success = request_slot(slot);
                    if (success) {
                        write_slot(slot, working_kmer);
                    }
                } while (!success && probe < local_size());

                    return false;     

                    

            }
    }   
        
    }
    */
return true;
}

///*
upcxx::future<> HashMap::send_buffer(upcxx::global_ptr<std::vector<kmer_pair>> remote_dst, int sending_rank, std::vector<kmer_pair> buf) {
  
  return upcxx::rpc(remote_dst.where(),
    [](const upcxx::global_ptr<std::vector<kmer_pair>> &dst, const int & sender_rank, upcxx::view<std::vector<kmer_pair>> &buf_in_rpc) {
      
      //
        std::vector<kmer_pair> * dst_recv_buff = dst.local(); 
        std::cout << sender_rank << std::endl;

        //for (auto kmer_it: buf_in_rpc) {
            //dst_recv_buff[sender_rank].push_back(*&kmer_it);
            //td::cout << kmer_it[0].kmer.get() << std::endl;
         //   std::cout << "a" << std::endl;

        //}
     
    },
    remote_dst, sending_rank, upcxx::make_view(buf.begin(), buf.end()));
}
//*/

/*

upcxx::future<> HashMap::send_buffer(upcxx::global_ptr<std::vector<kmer_pair>> receiving_proc, std::vector<kmer_pair> buf, int sending_proc) { 
                           
    return upcxx::rpc(receiving_proc,
                      // lambda to insert the key-value pair
                      [&](upcxx::global_ptr<std::vector<kmer_pair>>& dst, 
                                int sending_proc, std::vector<kmer_pair> buf) {
                        // insert into the local map at the target

                        std::vector<kmer_pair> * dst_recv_buff = receiving_proc.local(); 

                         for (auto& kmer : buf) {

                            dst_recv_buff[sending_proc].push_back(kmer);

                         }

                      }, receiving_proc, sending_proc, buf);
}
*/

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


HashMap::~HashMap() {
    // Need to destroy the atomic domain when the hashmap is destroyed
    ad.destroy();
    delete [] how_many_in_buffer;
}

