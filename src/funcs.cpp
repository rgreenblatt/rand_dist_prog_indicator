#include "funcs.h"

#include <chrono>
#include <thread>
#include <random>
#include <memory>
#include <functional>
#include <limits.h>
#include <assert.h>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/bind.hpp>
#include <iostream>
#include <tanhsinh.h>

void sleep_test(std::shared_ptr<std::atomic<unsigned int>> count, std::shared_ptr<std::atomic<double>> total, unsigned int dist_min, unsigned int dist_max, unsigned int total_sleep_count) {
    std::random_device rd;
    std::default_random_engine engine{rd()};
    std::uniform_int_distribution<unsigned int> dist(dist_min, dist_max);

    for(unsigned int i = 0; i < total_sleep_count; i++) {
        count->store(total_sleep_count - i);

        int sleep = dist(engine); 

        std::this_thread::sleep_for(std::chrono::milliseconds(sleep));

        total->store(total->load() + sleep);
    }
    count->store(0);
}

double boost_func_to_real(double x, const void * func) {
    
    auto boost_func_ptr = (boost::function<double(double)> * ) func;

    return (* boost_func_ptr)(x);
}


double num_integrate(boost::function<double(double)> func, double min, double max, double err_bound) {
    
    void * func_ptr = &func;

    return tanhsinh_quad(boost_func_to_real, func_ptr, min, max, 0.1, NULL, NULL); 
}

double cumlat_dist(double x, unsigned int remaining_samples) {
    double sum = 0;

    for(unsigned int k = 0; k <= remaining_samples; k++) {
        sum += std::pow(-1.0, k) * boost::math::binomial_coefficient<double>(remaining_samples, k) * std::pow(x - k, remaining_samples - 1) * boost::math::sign(x - k);
    } 

    return sum / (2.0 * boost::math::factorial<double>(remaining_samples - 1));
}

double shifted_cumlat_dist(double x, unsigned int remaining_samples, double min, double max, double adder) {
    return cumlat_dist((x - remaining_samples * min) / ((double) (max - min)), remaining_samples) / ((double) (max - min)) + adder;
}

double integrated_func(double x, std::vector<boost::function<double(double)>> shifted_cumlat_dist_filled_funcs, double err_bound) {
        
    double out = 1;
    double sum = 0;

    for(auto &shifted_cumlat_dist_filled : shifted_cumlat_dist_filled_funcs) {
        double f_at_x = 1.0 - num_integrate(shifted_cumlat_dist_filled, 0, x, err_bound);
        double df_at_x = -shifted_cumlat_dist_filled(x);
        out *= f_at_x;

        sum += df_at_x / f_at_x;
    }
    
    return x * -1.0 * out * sum;
}

double expected_val(std::vector<unsigned int> remaining_samples_per, std::vector<double> adders, double min, double max, double err_bound, unsigned int total_samples) {
    assert(remaining_samples_per.size() == adders.size());

    std::vector<boost::function<double(double)>> shifted_cumlat_dist_filled_funcs;

    for(size_t i = 0; i < adders.size(); i++) {
        shifted_cumlat_dist_filled_funcs.push_back(boost::bind(&shifted_cumlat_dist, ::_1, remaining_samples_per[i], min, max, adders[i]));
    } 
    
    auto integrated_filled = boost::bind(&integrated_func, ::_1, shifted_cumlat_dist_filled_funcs, err_bound);
    
    return num_integrate(integrated_filled, min * total_samples, max * total_samples, err_bound);
}

std::pair<unsigned int, unsigned int> get_min_max(std::vector<unsigned int> &vec) {
    std::pair<unsigned int, unsigned int> min_max(UINT_MAX, 0);
    for(auto &iter : vec) {
        if(iter < min_max.first) {
            min_max.first = iter;
        }
        if(iter > min_max.second) {
            min_max.second = iter;
        }
    }

    return min_max;
}
