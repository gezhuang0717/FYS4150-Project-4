#pragma once

#include <map>
#include <functional>
#include <vector>
#include <fstream>

namespace stat_utils{
    /**
     * @brief Write distribution to file
     * 
     * @param p Distribution
     * @param outfile File to write to
     */
    void write_distribution(std::map<double, double> p, std::ofstream &outfile);
    /**
     * @brief Calculate expected value from distribution
     * 
     * @param p Distribution
     * @return double 
     */
    double expected_value(std::map<double, double> p);
    /**
     * @brief Caluclate expected value from distrubution, under some transform
     * 
     * @param p Distribution
     * @param f Transform
     * @return double 
     */
    double expected_value(std::map<double, double> p, std::function<double(double)> f);
    /**
     * @brief Calculate expected value from N samples
     * 
     * @param samples Samples
     * @param N Samples to consider (N < size of samples)
     * @return double 
     */
    double expected_value(std::vector<int> samples, int N);
    /**
     * @brief Calculate expected value from N samples, under some transform
     * 
     * @param samples Samples
     * @param N Samples to consider (N < size of samples)
     * @param f Transform
     * @return double 
     */
    double expected_value(std::vector<int> samples, int N, std::function<double(int)> f);
    /**
     * @brief Find the distribution from N samples
     * 
     * @param samples Samples
     * @param N Samples to consider (N < size of samples)
     * @return std::map<double, double> 
     */
    std::map<double, double> distribution(std::vector<int> samples, int N);
        /**
     * @brief Find the distribution from N samples, under some transform
     * 
     * @param samples Samples
     * @param N Samples to consider (N < size of samples)
     * @param f Transform
     * @return std::map<double, double> 
     */
    std::map<double, double> distribution(std::vector<int> samples, int N, std::function<double(int)> f);
}

