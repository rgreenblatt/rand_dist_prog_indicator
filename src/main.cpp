#include "funcs.h"

#include <iostream>
#include <thread>
#include <chrono>
#include <memory>

int main() {

    std::vector<std::thread> sleep_threads;
   
    std::vector<std::shared_ptr<std::atomic<unsigned int>>> sleep_counts;
    std::vector<std::shared_ptr<std::atomic<double>>> total_sleep;

    unsigned int dist_min = 200;
    
    unsigned int dist_max = 650;
 
    unsigned int total_sleep_count = 10;

    unsigned int num_instances = 15;
    
    for(unsigned int i = 0; i < num_instances; i++) {
        sleep_counts.push_back(std::make_shared<std::atomic<unsigned int>>(total_sleep_count));
        total_sleep.push_back(std::make_shared<std::atomic<double>>(0));
        std::thread test(sleep_test, sleep_counts[i], total_sleep[i], dist_min, dist_max, total_sleep_count);
        sleep_threads.push_back(std::move(test)); 
    } 

    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout << std::endl;

    while(true) {
        std::vector<unsigned int> loaded_counts;
        std::vector<double> loaded_sleeps;

        for(auto &iter : sleep_counts) {
            loaded_counts.push_back(iter->load());
        }

        for(auto &iter : total_sleep) {
            loaded_sleeps.push_back(/*iter->load()*/ 0);
        }

        std::pair<unsigned int, unsigned int> min_max = get_min_max(loaded_counts);
    
        std::cout << /*"\x1b[A" << "\x1b[2K" <<*/ "Min count: " << min_max.first << " max count: " << min_max.second << " ETA: " << expected_val(loaded_counts, loaded_sleeps, dist_min, dist_max, 1, total_sleep_count) / 1000.0 - std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count() << std::endl; 

        std::this_thread::sleep_for(std::chrono::milliseconds(40));

        if(min_max.second == 0) {
            for(auto &iter : sleep_threads) {
                iter.join();
            }
            return 0; 
        }
    }

    return 1;
}
