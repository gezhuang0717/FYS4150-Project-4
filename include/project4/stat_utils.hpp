#pragma once

#include <map>
#include <functional>
#include <vector>
#include <fstream>

namespace stat_utils{
    void write_distribution(std::map<double, double> p, std::ofstream &outfile);
    double expected_value(std::map<double, double> p);
    double expected_value(std::map<double, double> p, std::function<double(double)> f);
    double expected_value(std::vector<int> samples, int N);
    double expected_value(std::vector<int> samples, int N, std::function<double(int)> f);
    std::map<double, double> distribution(std::vector<int> samples, int N);
    std::map<double, double> distribution(std::vector<int> samples, int N, std::function<double(int)> f);
}

