#include <gtest/gtest.h>
#include "funcs.h"
#include <chrono>
#include <random>
#include <boost/bind.hpp>

TEST(sleep, all) {
    std::default_random_engine engine{static_cast<unsigned int>(testing::UnitTest::GetInstance()->random_seed())};

    std::uniform_int_distribution<unsigned int> count_dist(3, 6);

    unsigned int total_sleep_count = count_dist(engine);

    std::shared_ptr<std::atomic<unsigned int>> count = std::make_shared<std::atomic<unsigned int>>(total_sleep_count); 

    std::shared_ptr<std::atomic<double>> total = std::make_shared<std::atomic<double>>(0.0); 

    std::uniform_int_distribution<unsigned int> sleep_time_dist(10, 15);

    unsigned int sleep_time = sleep_time_dist(engine);

    auto start_time = std::chrono::high_resolution_clock::now();
    sleep_test(count, total, sleep_time, sleep_time, total_sleep_count);
    auto end_time = std::chrono::high_resolution_clock::now();

    ASSERT_EQ(count->load(), 0);
    ASSERT_EQ(total->load(), total_sleep_count * sleep_time);
    ASSERT_NEAR(std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count(), total_sleep_count * sleep_time, 2);
}

double poly_func(double x, double a, double b, double c, double d, double e, double f) {
    return a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f;
}

TEST(integrate, value) {
    std::default_random_engine engine{static_cast<unsigned int>(testing::UnitTest::GetInstance()->random_seed())};
    
    std::uniform_real_distribution<double> coefs_dist(0.0, 100.0);
    
    double a = 0;
    double b = coefs_dist(engine);
    double c = coefs_dist(engine);
    double d = coefs_dist(engine);
    double e = coefs_dist(engine);
    double f = coefs_dist(engine);

    std::uniform_real_distribution<double> min_dist(0.0, 100.0);
    std::uniform_real_distribution<double> max_dist(110.0, 200.0);

    double min = min_dist(engine);
    double max = max_dist(engine);

    auto integrated = boost::bind(&poly_func, ::_1, b / 5.0, c / 4.0, d / 3.0, e / 2.0, f, 0.0);

    double actual = integrated(max) - integrated(min);

    ASSERT_NEAR(num_integrate(boost::bind(&poly_func, ::_1, a, b, c, d, e, f), min, max, 0.1), actual, 0.1);

}

struct cumlat_dist_tests {
    double expected = 0.0;
    double x = 0.0;
    unsigned int remaining_samples = 0;
};

TEST(cumlat_dist, value) {
    std::vector<cumlat_dist_tests> tests = {
    {0.5, 1.5, 2},
    {0.125, 0.5, 3},
    {0.361503940661, 5.1, 9},
    {0.545024166667, 3.1, 6},
    };
    for(auto &test : tests) {
        ASSERT_NEAR(cumlat_dist(test.x, test.remaining_samples), test.expected, test.expected * 0.001);
    }
}

struct shifted_cumlat_dist_tests {
    double expected = 0.0;
    double x = 0.0;
    unsigned int remaining_samples = 0;
    double min = 0.0;
    double max = 0.0;
};

TEST(shifted_cumlat_dist, value) {
    std::vector<shifted_cumlat_dist_tests> tests = {
    {0.5, 2.5, 2, 1.0, 2.0},
    {0.108201505412, 17.65, 3, 3.0, 9.6},
    {0.0278264853498, 53.65, 9, 5, 8},
    };

    for(auto &test : tests) {
        ASSERT_NEAR(shifted_cumlat_dist(test.x, test.remaining_samples, test.min, test.max), test.expected, test.expected * 0.001);
    }
}

struct expected_val_tests {
    double expected = 0.0;
    std::vector<unsigned int> remaining_samples_per;
    std::vector<double> adders;
    double min = 0.0;
    double max = 0.0;
    unsigned int steps = 1000;
    unsigned int total_samples = 100;
};

TEST(expected_val, value) {
    std::vector<expected_val_tests> tests = {
    {57.0300758469, {9, 9}, {0, 0}, 5, 8, 1000, 9},
    };

    for(auto &test : tests) {
        ASSERT_NEAR(expected_val(test.remaining_samples_per, test.adders, test.min, test.max, test.steps, test.total_samples), test.expected, test.expected * 0.05);
    }
}


int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
