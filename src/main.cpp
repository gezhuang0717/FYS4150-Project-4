#include "project4/ising_model.hpp"
#include "project4/stat_utils.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include <memory>
#include <string>
#include <chrono>

using namespace std;

void sample(vector<int> &sampled_energy, vector<int> &sampled_magnetization_abs, const int iters, int L, double T, int seed, int burn_in_time = 1000, bool random_spins = true){
    IsingModel model(L, T, random_spins, seed);
    for (int i = -burn_in_time; i < iters; i++){
        model.metropolis();
        if (i < 0) continue;
        sampled_energy.push_back(model.get_energy());
        sampled_magnetization_abs.push_back(abs(model.get_magnetization()));
    }
}


/**
 * @brief Produces sampled distribution of &epsilon; and |m|
 * 
 * @param N Numbers of MCMC-iterations
 * @param L Size of the Lattice (will have L * L elements)
 * @param T temperature
 * @param burn_in_time Number of iterations to discard when producing estimates 
 */
void write_distributions(const int iters, const int L, double T, int seed, string filename_epsilon, string filename_m_abs, int burn_in_time = 1000){ // note we don't have any indications on what the burn_in_time is yet
    // IsingModel model(L, T, true, seed); // unused?
    const int N = L * L;
    vector<int> sampled_energy;
    vector<int> sampled_magnetization_abs;
    sample(sampled_energy, sampled_magnetization_abs, iters, L, T, seed);
    auto scale = [N](int x){return (double)x / N;};
    map<double, double> buckets_epsilon = stat_utils::distribution(sampled_energy, iters, scale);
    map<double, double> buckets_m_abs = stat_utils::distribution(sampled_magnetization_abs, iters, scale);
    ofstream outfile_epsilon(filename_epsilon);
    outfile_epsilon << "epsilon,p" << endl;
    stat_utils::write_distribution(buckets_epsilon, outfile_epsilon);
    ofstream outfile_m_abs(filename_m_abs);
    outfile_epsilon << "|m|,p" << endl;
    stat_utils::write_distribution(buckets_m_abs, outfile_m_abs);
}


void values(vector<int> sampled_energy, vector<int> sampled_magnetization_abs, int sample_size, int L, double T, 
double &expected_epsilon, double &expected_m_abs, double &c_v, double &chi){
    int N = L * L;
    auto scale = [N](int x){return (double)x / N;};
    auto square = [](int x){return x * x;};
    expected_epsilon = stat_utils::expected_value(sampled_energy, sample_size, scale);
    expected_m_abs = stat_utils::expected_value(sampled_magnetization_abs, sample_size, scale);
    double expected_energy = stat_utils::expected_value(sampled_energy, sample_size);
    double expected_energy_sq = stat_utils::expected_value(sampled_energy, sample_size, square);
    double expected_magnetization_abs = stat_utils::expected_value(sampled_magnetization_abs, sample_size);
    double expected_magnetization_sq = stat_utils::expected_value(sampled_magnetization_abs, sample_size, square);
    c_v = (1. / N) * (1. / (T * T)) * (expected_energy_sq - expected_energy * expected_energy);
    chi = (1. / N) * (1. / T) * (expected_magnetization_sq - expected_magnetization_abs * expected_magnetization_abs);
}



/**
 * @brief estimates <&epsilon;>, <|m|>, C_v and &chiand writes them to the outfile in that order, all after the temperature
 * @param L problem size
 * @param T temperature
 * @param outfile csv-file to write results
 */
void write_values_to_file(int L, double T, int seed, ofstream &outfile){
    int sample_size = 10000;
    vector<int> sampled_energy;
    vector<int> sampled_magnetization_abs;
    sample(sampled_energy, sampled_magnetization_abs, sample_size, L, T, seed, 1000);
    double expected_epsilon, expected_m_abs, c_v, chi;
    values(sampled_energy, sampled_magnetization_abs, sample_size, L, T, expected_epsilon, expected_m_abs, c_v, chi);
    outfile << T << "," << expected_epsilon << "," << expected_m_abs << "," << c_v << "," << chi << endl;
}

void find_burn_in_time(int N, int L, double T, int seed, bool random_spins=true){
    vector<int> sampled_E;
    vector<int> sampled_M_abs;
    vector<double> E_avg;
    vector<double> M_abs_avg;
    vector<int> N_vector;
    sample(sampled_E, sampled_M_abs, N, L, T, seed, 0, random_spins);
    for (int i=0; i<N; i++){
        double expected_eps, expected_m_abs, c_v, chi;
        values(sampled_E, sampled_M_abs, i+1, L, T, expected_eps, expected_m_abs, c_v, chi);
        E_avg.push_back(expected_eps);
        M_abs_avg.push_back(expected_m_abs);
        N_vector.push_back(i+1);
    }

    ofstream burn_in_csv;
    string spins_are_random;
    if (random_spins){
        spins_are_random = "random";
    }
    else{
        spins_are_random = "nonrandom";
    }
    string filename = "output/burn_in_T_" + to_string(T) + "_" + spins_are_random + ".csv";
    burn_in_csv.open(filename);
    burn_in_csv << "N,expected_E,expected_M\n";

    for (int i=0; i<N; i++){
        burn_in_csv << N_vector[i] << "," << E_avg[i] 
                    << "," << M_abs_avg[i] << "\n";
    }
    burn_in_csv.close();
}

/**
 * @brief Testing for convergence against analytical results in the 2x2 case
 * 
 * @return int number of samples needed for convergence
 */
int test2x2(){
    int seed = 3875623;
    const double tol = 1e-2;
    vector<int> sampled_energy;
    vector<int> sampled_magnetization_abs;
    double T = 1.;
    int max_sample_size = 10000;
    sample(sampled_energy, sampled_magnetization_abs, max_sample_size, 2, T, seed, 0);
    double expected_epsilon, expected_m_abs, c_v, chi;
    double analytical_expected_epsilon = -1.99598208593669, analytical_expected_m_abs = 0.9986607327485997, analytical_c_v = 0.032082331864287994, analytical_chi = 0.004010739516227435;
    int using_sample_size = 1;
    do {
        values(sampled_energy, sampled_magnetization_abs, using_sample_size, 2, T, expected_epsilon, expected_m_abs, c_v, chi);
        using_sample_size++;
    } while((abs(expected_epsilon - analytical_expected_epsilon) > tol 
            or abs(expected_m_abs - analytical_expected_m_abs) > tol 
            or abs(c_v - analytical_c_v) > tol 
            or abs(chi - analytical_chi) > tol)
            and using_sample_size < max_sample_size);
    return using_sample_size;
}

void timing_parallel_vs_serial(int L, double T) {
    int repeats = 10;

    double T_min = 2.1;
    double T_max = 2.1;
    int steps = 48;
    double dT = (T_max - T_min) / steps;
    int seed = 42;

    double total_time_serial = 0;
    double total_time_parallel = 0;

    cout << "Running parallel: " << flush;
    for (int i = 0; i < repeats; i++) {
        auto start = chrono::high_resolution_clock::now();
        ofstream outfile1("output/values_L=" + to_string(L) + ".csv");
        #pragma omp parallel for
        for (int i = 0; i < steps; i++){
            double T = T_min + i * dT;
            write_values_to_file(L, T, seed, outfile1);
        }
        outfile1.close();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_parallel = end - start;
        total_time_parallel += diff_parallel.count();
        cout << "." << flush;
    }


    cout << "Running serial: " << flush;
    for (int i = 0; i < repeats; i++) {
        auto start = chrono::high_resolution_clock::now();
        ofstream outfile2("output/values_L=" + to_string(L) + ".csv");
        for (int i = 0; i < steps; i++){
            double T = T_min + i * dT;
            write_values_to_file(L, T, seed, outfile2);
        }
        outfile2.close();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_serial = end - start;
        total_time_serial += diff_serial.count();
        cout << "." << flush;
    }

    ofstream timingfile("output/timing.txt");
    timingfile << "Average timing for L=" << L << ", T=" << T << " and steps=" << steps << ", with " << repeats << " repeats" << endl;

    timingfile << "Average time for parallel: " << total_time_parallel / repeats << endl;
    timingfile << "Average time for serial: " << total_time_serial / repeats << endl;
    timingfile << "Parallel is " << total_time_serial / total_time_parallel << " times faster than serial" << endl;

    timingfile.close();
}

int main(){
    timing_parallel_vs_serial(40, 2);
    return 0;

    cout << "Testing for convergence against analytical results in the 2x2 case. Needed sample size: " << test2x2() << endl;
    int seed = 456788;
    double T_min = 2.1;
    int steps = 48;

    for (int L = 20; L <= 100; L += 20) {
        cout << "Testing for " << L << "x" << L << endl;
        double dT = (2.4 - T_min) / steps;
        ofstream outfile("output/values_L=" + to_string(L) + ".csv");
        #pragma omp parallel for
        for (int i = 0; i < steps; i++){
            double T = T_min + i * dT;
            write_values_to_file(L, T, seed++, outfile);
        }
        outfile.close();
    }

    
    find_burn_in_time(2000, 20, 1, seed++, false);
    find_burn_in_time(1000, 20, 1, seed++, true);
    find_burn_in_time(6000, 20, 2.4, seed++, false);
    find_burn_in_time(6000, 20, 2.4, seed++, true);

    return 0;
}
