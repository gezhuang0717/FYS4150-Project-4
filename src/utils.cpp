#include <iostream>
#include <map>
#include "project4/stat_utils.hpp"

double I(int x) {return x;}

void stat_utils::print_map(std::map<int, double> p){
    for (auto const& bucket : p){
        std::cout << bucket.first << ": " << bucket.second << std::endl;
    }
}


double stat_utils::expected_value(std::map<int, double> p, std::function<double(int)> f){
    double expected_value = 0.;
    for (auto const& bucket : p){
        expected_value += f(bucket.first) * bucket.second;
    }
    return expected_value;
}

double stat_utils::expected_value(std::map<int, double> p){
    return expected_value(p, I);
}

double stat_utils::expected_value(int *samples, int N, std::function<double(int)> f){
    double sum = 0;
    for (int i = 0; i < N; i++){
        sum += f(samples[i]);
    }
    return sum / N;
}

double stat_utils::expected_value(int *samples, int N){
    return expected_value(samples, N, I);
}

std::map<int, double> stat_utils::distribution(int *samples, int N, std::function<double(int)> f){
    std::map<int, double> buckets;
    for (int i = 0; i < N; i++){
        float sample = f(samples[i]);
        buckets.emplace(sample, 0);
        buckets[sample] += 1./N;
    }
    return buckets;
}

std::map<int, double> stat_utils::distribution(int *samples, int N){
    return distribution(samples, N, I);
}
