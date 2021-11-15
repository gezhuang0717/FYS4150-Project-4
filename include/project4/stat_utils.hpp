#pragma once

#include <map>
#include <functional>

namespace stat_utils{
    void print_map(std::map<int, float> p);
    float expected_value(std::map<int, float> p);
    float expected_value(std::map<int, float> p, std::function<float(int)> f);
    float expected_value(int* samples, int N);
    float expected_value(int* samples, int N, std::function<float(int)> f);
    std::map<int, float> distribution(int* samples, int N);
    std::map<int, float> distribution(int* samples, int N, std::function<float(int)> f);
}

