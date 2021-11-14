#include <iostream>
#include <map>
#include "project4/stat_utils.hpp"

int I(int x) {return x;}

void stat_utils::print_map(std::map<int, float> p){
    for (auto const& bucket : p){
        std::cout << bucket.first << ": " << bucket.second << std::endl;
    }
}


float stat_utils::expected_value(std::map<int, float> p, std::function<int(int)> f){
    float expected_value = 0.;
    for (auto const& bucket : p){
        expected_value += f(bucket.first) * bucket.second;
    }
    return expected_value;
}

float stat_utils::expected_value(std::map<int, float> p){
    return expected_value(p, I);
}

float stat_utils::expected_value(int *samples, int N, std::function<int(int)> f){
    int sum = 0;
    for (int i = 0; i < N; i++){
        sum += f(samples[i]);
    }
    return (float)sum / N;
}

float stat_utils::expected_value(int *samples, int N){
    return expected_value(samples, N, I);
}

std::map<int, float> stat_utils::distribution(int *samples, int N, std::function<int(int)> f){
    std::map<int, float> buckets;
    for (int i = 0; i < N; i++){
        int sample = f(samples[i]);
        buckets.emplace(sample, 0);
        buckets[sample] += 1./N;
    }
    return buckets;
}

std::map<int, float> stat_utils::distribution(int *samples, int N){
    return distribution(samples, N, I);
}
