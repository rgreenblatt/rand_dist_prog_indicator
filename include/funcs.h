#include <atomic>
#include <vector>
#include <memory>
#include <boost/function.hpp>

void sleep_test(std::shared_ptr<std::atomic<unsigned int>> count, std::shared_ptr<std::atomic<double>> total, unsigned int dist_min, unsigned int dist_max, unsigned int total_sleep_count);

double num_integrate(boost::function<double(double)> func, double min, double max, double err_bound);

double cumlat_dist(double x, unsigned int remaining_samples);

double shifted_cumlat_dist(double x, unsigned int remaining_samples, double min, double max, double adder);

double integrated_func(double x, std::vector<boost::function<double(double)>> shifted_cumlat_dist_filled_funcs, double err_bound);

double expected_val(std::vector<unsigned int> remaining_samples_per, std::vector<double> adders, double min, double max, double err_bound, unsigned int total_samples);

std::pair<unsigned int, unsigned int> get_min_max(std::vector<unsigned int> &vec);
