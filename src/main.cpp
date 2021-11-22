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
#include <sstream>

using namespace std;

bool has_flag(const std::string& option, char** begin, char** end){
    return std::find(begin, end, option) != end;
}

void print_help_message() {
    cout << "Usage" << endl;
    cout << "\t./runner [flags]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "\t-h\tShow this help message" << endl;
    cout << "\t-t\tTest implementation" << endl;
    cout << "\t-b\tFinds burn-in time" << endl;
    cout << "\t-w\tWrites samples to file. Provide L, T, and seed" << endl;
    cout << "\t-s\tFinds values for L in 20, 40, 60, 80, 100, 120, 140, 160 for T in the range [2.1, 2.4]" << endl;
    cout << "\t-z\tZooms in and finds values. Provide L, T_min, T_max and seed as system argunments" << endl;
}

/*
 * @brief Samples from IsingModel, once per MC-cycle
 * 
 * @param sampled_energy Destination of the sampled energy
 * @param sampled_magnetization_abs Destination of the samped absouloute magnetization
 * @param iters Number of samples to make
 * @param L Size of the IsingModel
 * @param T Temperature of the IsingModel
 * @param seed Seed for RNG
 * @param burn_in_time Number of iterations to discard
 * @param random_spins If True, starting spin state of the IsingModel is chosen randomly. Else all spins are poining in the same direction
 */
void sample(vector<int> &sampled_energy, vector<int> &sampled_magnetization_abs, const int iters, int L, double T, int seed, int burn_in_time = 10000, bool random_spins = true){
    IsingModel model(L, T, random_spins, seed);
    for (int i = -burn_in_time; i < iters; i++){
        model.metropolis();
        if (i < 0) continue;
        sampled_energy.push_back(model.get_energy());
        sampled_magnetization_abs.push_back(abs(model.get_magnetization()));
    }
}


/**
 * @brief writes samples of &epsilon; and |m|
 * 
 * @param N Numbers of MC-iterations
 * @param L Size of the Lattice (will have L * L elements)
 * @param T temperature
 * @param burn_in_time Number of iterations to discard when producing estimates 
 */
void write_samples(const int iters, const int L, double T, int seed, int burn_in_time = 10000){
    const int N = L * L;
    vector<int> sampled_energy;
    vector<int> sampled_magnetization_abs;
    sample(sampled_energy, sampled_magnetization_abs, iters, L, T, seed);
    auto scale = [N](int x){return (double)x / N;};
    ostringstream out;
    out.precision(1);
    out << fixed << T;
    ofstream outfile("output/samples_L=" + to_string(L) + "_T=" + out.str() + ".csv");
    outfile << "epsilon,|m|" << endl;
    for (int i = 0; i < iters; i++){
        outfile << scale(sampled_energy[i]) << "," << sampled_magnetization_abs[i] << endl;
    }
    outfile.close();
}

/**
 * @brief From sampled values, estimates <&epsilon;> <|m|>, C_v and &chi;
 * 
 * @param sampled_energy Energy sampled from IsingModel
 * @param sampled_magnetization_abs Absuloute magnetization sampled from IsingModel
 * @param sample_size Size of the samples to consider. Note this could be less than the number of samples available, if we would like to estimate using only the first samples
 * @param L Size of the IsingModel
 * @param T Temperature of the IsingModel
 * @param expected_epsilon Destinatipn of <&epsilon;>
 * @param expected_m_abs Destination of <|m|>
 * @param c_v Destination of C_v
 * @param chi Destination of &chi;
 */
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
    int sample_size = 1000000;
    vector<int> sampled_energy;
    vector<int> sampled_magnetization_abs;
    sample(sampled_energy, sampled_magnetization_abs, sample_size, L, T, seed, 30000);
    double expected_epsilon, expected_m_abs, c_v, chi;
    values(sampled_energy, sampled_magnetization_abs, sample_size, L, T, expected_epsilon, expected_m_abs, c_v, chi);
    outfile << T << "," << expected_epsilon << "," << expected_m_abs << "," << c_v << "," << chi << endl;
}

/**
 * @brief Trying to estimate the burn in time in the IsingModel
 * 
 * @param N TODO: Rename?
 * @param L Size of the IsingModel
 * @param T Temperature of the IsingModel
 * @param seed seed for RNG
 * @param random_spins f True, starting spin state of the IsingModel is chosen randomly. Else all spins are poining in the same direction
 */
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
    string filename = "output/burn_in_L_" + to_string(L) + "_T_" 
                    + to_string(T) + "_" + spins_are_random + ".csv";
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
int test2x2(int seed){
    const double tol = 1e-3;
    vector<int> sampled_energy;
    vector<int> sampled_magnetization_abs;
    double T = 2.;
    int max_sample_size = 10000;
    sample(sampled_energy, sampled_magnetization_abs, max_sample_size, 2, T, seed, 0);
    double expected_epsilon, expected_m_abs, c_v, chi;
    double analytical_expected_epsilon = -1.8008253628497959, analytical_expected_m_abs = 0.9337091730054017, analytical_c_v = 0.3610959875477656, analytical_chi = 0.09079837108784634;
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

/**
 * @brief Performs timing to compare parallel and serial code
 * 
 * @param L Size of the IsingModel
 * @param T Temperature of the IsingModel
 */
void timing_parallel_vs_serial(int L, double T) {
    int repeats = 10;

    double T_min = 2.1;
    double T_max = 2.1;
    int steps = 48;
    double dT = (T_max - T_min) / steps;
    int seed = 42;

    double total_time_serial = 0;
    double total_time_parallel = 0;

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
    }


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
    }

    ofstream timingfile("output/timing.txt");
    timingfile << "Average timing for L=" << L << ", T=" << T << " and steps=" << steps << ", with " << repeats << " repeats" << endl;

    timingfile << "Average time for parallel: " << total_time_parallel / repeats << endl;
    timingfile << "Average time for serial: " << total_time_serial / repeats << endl;
    timingfile << "Parallel is " << total_time_serial / total_time_parallel << " times faster than serial" << endl;

    timingfile.close();
}

/**
 * @brief Writes estimated values from IsingModels within the specified range of temperatures
 * Important: If this is ran on multiple cores, you may not assume the file is sorted by temperature!
 * 
 * @param T_min The minimum temperature to estimate for
 * @param T_max The maximum temperature to estimate for (non inclusive)
 * @param L The size of the IsingModel
 * @param steps The number of temperatures between T_min and T_max
 * @param seed Seed for RNG
 * @param filename Filename of the outputfile
 */
void look_between_temperatures(double T_min, double T_max, int L, int steps, int &seed, string filename){
    cout << "Testing for " << L << "x" << L << endl;
    double dT = (T_max - T_min) / steps;
    ofstream outfile(filename);
    outfile << "T,<epsilon>,<|m|>,C_v,chi" << endl;
    #pragma omp parallel for
    for (int i = 0; i < steps; i++){
        double T = T_min + i * dT;
        write_values_to_file(L, T, seed++, outfile);
    }
    outfile.close();
}



int main(int argc, char *argv[]){
    if (has_flag("-h", argv, argv + argc)) print_help_message();
    else if (has_flag("-t", argv, argv + argc)) {
        int seed = 3875623;
        cout << "Testing for convergence against analytical results in the 2x2 case. Needed sample size: " << test2x2(seed) << endl;
        timing_parallel_vs_serial(20, 2.);
    }
    else if (has_flag("-b", argv, argv + argc)) { 
        int seed =  23344;
        // find burn-in time for L=20
        find_burn_in_time(40000, 20, 1, seed++, false);
        find_burn_in_time(40000, 20, 1, seed++, true);
        find_burn_in_time(40000, 20, 2.4, seed++, false);
        find_burn_in_time(40000, 20, 2.4, seed++, true);

        // find burn-in time for L=100
        find_burn_in_time(40000, 100, 1, seed++, false);
        find_burn_in_time(40000, 100, 1, seed++, true);
        find_burn_in_time(40000, 100, 2.4, seed++, false);
        find_burn_in_time(40000, 100, 2.4, seed++, true);

    }
    else if (has_flag("-w", argv, argv + argc)){
        if (argc < 4){
            cout << "Please include L, T and seed" << endl;
            return 1;
        }
        int L = atoi(argv[2]);
        double T = atof(argv[3]);
        int seed =atoi(argv[4]);
        write_samples(100000, L, T, seed);
    }
    
    else if (has_flag("-s", argv, argv + argc)){
        int seed = 456788;
        int steps = 24;
        double T_min = 2.1;
        double T_max = 2.4;
        for (int L = 20; L <= 140; L += 20) {
            look_between_temperatures(T_min, T_max, L, steps, seed, "output/values_L=" + to_string(L) + ".csv");
        }
    }
    else if (has_flag("-z", argv, argv + argc)){
        int steps = 24;
        if (argc < 5){
             cout << "Please include L, T_min, T_max and seed" << endl;
            return 1;
        }
        int L = atoi(argv[2]);
        double T_min = atof(argv[3]);
        double T_max = atof(argv[4]);
        int seed = atoi(argv[5]);
        look_between_temperatures(T_min, T_max, L, steps, seed, "output/values_zoom_L=" + to_string(L) + ".csv");
    }
    return 0;
}
