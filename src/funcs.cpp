#include "funcs.h"

#include <chrono>
#include <thread>
#include <random>
#include <memory>
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

double num_integrate(boost::function<double(double)> func, double min, double max, double err_bound, bool &infinite_error_interval) {
    void * func_ptr = &func;
   
    double err_out = 0;
    unsigned count = 0;
 
    double retval = tanhsinh_quad(boost_func_to_real, func_ptr, min, max, err_bound, &err_out, &count); 

    infinite_error_interval = false;

    if(boost::math::isinf(err_out)) {
        retval = 0;
        infinite_error_interval = true;
    }

    return retval; 
}

double cumlat_dist(double x, unsigned int remaining_samples) {
    double sum = 0;
    
    for(unsigned int k = 0; k <= remaining_samples; k++) {
        sum += std::pow(-1.0, k) * boost::math::binomial_coefficient<double>(remaining_samples, k) * std::pow(x - k, remaining_samples - 1) * boost::math::sign(x - k);
    } 
    
    return sum / (2.0 * boost::math::factorial<double>(remaining_samples - 1));
}

double shifted_cumlat_dist(double x, unsigned int remaining_samples, double min, double max, double adder) {
    return cumlat_dist((x - remaining_samples * min - adder) / ((double) (max - min)), remaining_samples) / ((double) (max - min));
}

double integrated_func(double x, std::vector<boost::function<double(double)>> shifted_cumlat_dist_filled_funcs, double err_bound, bool min_or_max, std::vector<double> individual_maxes, std::vector<double> individual_mins) {
    assert(shifted_cumlat_dist_filled_funcs.size() ==  individual_maxes.size() && individual_maxes.size() == individual_mins.size());

    double out = 1;
    double sum = 0;

    for(size_t i = 0; i < shifted_cumlat_dist_filled_funcs.size(); i++) {
        double f_at_x = 1;
        double df_at_x = 0;

        double individual_x = x > individual_maxes[i] ? individual_maxes[i] : x;
        
        if(min_or_max) {
            
            bool err;

            f_at_x = 1.0 - num_integrate(shifted_cumlat_dist_filled_funcs[i], individual_mins[i], individual_x, err_bound / 10.0, err);

            if(err) {
                continue; 
            }
            
            df_at_x = -shifted_cumlat_dist_filled_funcs[i](individual_x);
        }
        else {
            bool err;

            f_at_x = num_integrate(shifted_cumlat_dist_filled_funcs[i], individual_mins[i], individual_x, err_bound / 10.0, err);

            if(err) {
                continue; 
            }
            
            df_at_x = shifted_cumlat_dist_filled_funcs[i](individual_x);
        }
        if(f_at_x != 0) {
            out *= f_at_x;

            sum += df_at_x / f_at_x;
        }
    }

    return x * out * sum * (min_or_max ? -1.0 : 1.0);
}

double expected_val(std::vector<unsigned int> remaining_samples_per, std::vector<double> adders, double min, double max, double err_bound, bool min_or_max) {
    assert(remaining_samples_per.size() == adders.size());
    assert(adders.size() > 0);

    std::vector<boost::function<double(double)>> shifted_cumlat_dist_filled_funcs;

    double minimum_min = DBL_MAX;
    double maximum_max = DBL_MIN;

    bool no_remaining_samples = false;
    double retval = min_or_max ? DBL_MAX : DBL_MIN;

    std::vector<double> individual_maxes;
    std::vector<double> individual_mins;

    for(size_t i = 0; i < adders.size(); i++) {
        if(remaining_samples_per[i] == 0) {
            no_remaining_samples = true;
            if(min_or_max) {
                retval = adders[i] < retval ? adders[i] : retval;
            }
            else {
                retval = adders[i] > retval ? adders[i] : retval;
            }
        }
        else if(!min_or_max) {
            no_remaining_samples = false;
        }

        if(!no_remaining_samples) {
            double individual_min = remaining_samples_per[i] * min + adders[i];
            double individual_max = remaining_samples_per[i] * max + adders[i];

            maximum_max = individual_max > maximum_max ? individual_max : maximum_max;
            minimum_min = individual_min < minimum_min ? individual_min : minimum_min;

            shifted_cumlat_dist_filled_funcs.push_back(boost::bind(&shifted_cumlat_dist, ::_1, remaining_samples_per[i], min, max, adders[i]));

            individual_mins.push_back(individual_min);
            individual_maxes.push_back(individual_max);
        }
    } 

    if((!no_remaining_samples && min_or_max) || ((shifted_cumlat_dist_filled_funcs.size() > 0) && !min_or_max)) {
        auto integrated_filled = boost::bind(&integrated_func, ::_1, shifted_cumlat_dist_filled_funcs, err_bound, min_or_max, individual_maxes, individual_mins);
       
        bool err;
 
        retval = num_integrate(integrated_filled, minimum_min, maximum_max, err_bound, err);

        assert(!err);
    }
    
    return retval;
}

std::pair<double, double> get_min_max(std::vector<double> &vec) {
    std::pair<double, double> min_max(DBL_MAX, DBL_MIN);
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
