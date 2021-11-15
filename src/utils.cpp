#include <iostream>
#include <map>
#include "project4/stat_utils.hpp"

float I(int x) {return x;}

void stat_utils::print_map(std::map<int, double> p){
    for (auto const& bucket : p){
        std::cout << bucket.first << ": " << bucket.second << std::endl;
    }
}


<<<<<<< HEAD
float stat_utils::expected_value(std::map<int, float> p, std::function<float(int)> f){
    float expected_value = 0.;
=======
double stat_utils::expected_value(std::map<int, double> p, std::function<int(int)> f){
    double expected_value = 0.;
>>>>>>> 53e76a0c940bec23214d6f3a9b1d7445615cc574
    for (auto const& bucket : p){
        expected_value += f(bucket.first) * bucket.second;
    }
    return expected_value;
}

double stat_utils::expected_value(std::map<int, double> p){
    return expected_value(p, I);
}

<<<<<<< HEAD
float stat_utils::expected_value(int *samples, int N, std::function<float(int)> f){
    float sum = 0;
    for (int i = 0; i < N; i++){
        sum += f(samples[i]);
    }
    return sum / N;
=======
double stat_utils::expected_value(int *samples, int N, std::function<int(int)> f){
    int sum = 0;
    for (int i = 0; i < N; i++){
        sum += f(samples[i]);
    }
    return (double)sum / N;
>>>>>>> 53e76a0c940bec23214d6f3a9b1d7445615cc574
}

double stat_utils::expected_value(int *samples, int N){
    return expected_value(samples, N, I);
}

<<<<<<< HEAD
std::map<int, float> stat_utils::distribution(int *samples, int N, std::function<float(int)> f){
    std::map<int, float> buckets;
=======
std::map<int, double> stat_utils::distribution(int *samples, int N, std::function<int(int)> f){
    std::map<int, double> buckets;
>>>>>>> 53e76a0c940bec23214d6f3a9b1d7445615cc574
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
