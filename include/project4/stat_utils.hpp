#pragma once

#include <map>
#include <functional>

namespace stat_utils{
    void print_map(std::map<double, double> p);
    double expected_value(std::map<double, double> p);
    double expected_value(std::map<double, double> p, std::function<double(double)> f);
    double expected_value(int* samples, int N);
    double expected_value(int* samples, int N, std::function<double(int)> f);
    std::map<double, double> distribution(int* samples, int N);
    std::map<double, double> distribution(int* samples, int N, std::function<double(int)> f);
}

