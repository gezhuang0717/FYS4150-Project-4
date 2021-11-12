#pragma once

#include <map>


namespace map_utils{
    void print_map(std::map<int, float> p);
    float expected_value(std::map<int, float> p);
    std::map<int, float> distribution(int* samples, int N);
}

