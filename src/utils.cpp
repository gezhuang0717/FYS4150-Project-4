#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "project4/stat_utils.hpp"

double I(float x){return x;}

void stat_utils::write_distribution(std::map<double, double> p, std::ofstream &outfile){
    for (auto const& bucket : p){
        outfile << bucket.first << ": " << bucket.second << std::endl;
    }
    outfile.close();
}


double stat_utils::expected_value(std::map<double, double> p, std::function<double(double)> f){
    double expected_value = 0.;
    for (auto const& bucket : p){
        expected_value += f(bucket.first) * bucket.second;
    }
    return expected_value;
}

double stat_utils::expected_value(std::map<double, double> p){
    return expected_value(p, I);
}

double stat_utils::expected_value(std::vector<int> samples, int N, std::function<double(int)> f){
    double sum = 0;
    for (int i = 0; i < N; i++){
        sum += f(samples[i]);
    }
    return sum / N;
}

double stat_utils::expected_value(std::vector<int> samples, int N){
    return expected_value(samples, N, I);
}

std::map<double, double> stat_utils::distribution(std::vector<int> samples, int N, std::function<double(int)> f){
    std::map<double, double> buckets;
    for (int i = 0; i < N; i++){
        float sample = f(samples[i]);
        buckets.emplace(sample, 0);
        buckets[sample] += 1./N;
    }
    return buckets;
}

std::map<double, double> stat_utils::distribution(std::vector<int> samples, int N){
    return distribution(samples, N, I);
}