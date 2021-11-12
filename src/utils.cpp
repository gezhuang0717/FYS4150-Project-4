#include <iostream>
#include <map>
#include "project4/map_utils.hpp"

void map_utils::print_map(std::map<int, float> p){
    for (auto const& bucket : p){
        std::cout << bucket.first << ": " << bucket.second << std::endl;
    }
}

float map_utils::expected_value(std::map<int, float> p){
    float expected_value = 0.;
    for (auto const& bucket : p){
        expected_value += bucket.first * bucket.second;
    }
    return expected_value;
}

std::map<int, float> map_utils::distribution(int *samples, int N){
    std::map<int, float> buckets;
    for (int i = 0; i < N; i++){
        int sample = samples[i];
        buckets.emplace(sample, 0);
        buckets[sample] += 1./N;
    }
    return buckets;
}