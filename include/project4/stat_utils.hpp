#pragma once

#include <map>
#include <functional>

namespace stat_utils{
    void print_map(std::map<int, double> p);
    double expected_value(std::map<int, double> p);
    double expected_value(std::map<int, double> p, std::function<int(int)> f);
    double expected_value(int* samples, int N);
    double expected_value(int* samples, int N, std::function<int(int)> f);
    std::map<int, double> distribution(int* samples, int N);
    std::map<int, double> distribution(int* samples, int N, std::function<int(int)> f);
}

