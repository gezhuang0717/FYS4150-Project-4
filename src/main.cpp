#include "project4/ising_model.hpp"
#include "project4/stat_utils.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>

using namespace std;

void find_burn_in_time(int N, int L=20, double T=1){
    IsingModel model(L, T, true);
    
    int sampled_E[N];
    int sampled_M[N];
    for (int i=0; i<N; i++){
        sampled_E[i] = model.get_energy();
        sampled_M[i] = model.get_magnetization();
        model.metropolis();
    }
    vector<int> burn_in_time;
    vector<double> expected_E_burn_in;
    vector<double> expected_M_burn_in;
    for (int burn_in=0; burn_in<N; burn_in+=max(1, N/100)){
        map<int, double> buckets_E = stat_utils::distribution(sampled_E + burn_in, N - burn_in);
        map<int, double> buckets_M = stat_utils::distribution(sampled_M + burn_in, N - burn_in);
        expected_E_burn_in.push_back(stat_utils::expected_value(buckets_E));
        expected_M_burn_in.push_back(stat_utils::expected_value(buckets_M));
        burn_in_time.push_back(burn_in);
    }

    ofstream burn_in_csv;
    burn_in_csv.open("output/burn_in.csv");
    burn_in_csv << "burn_in,expected_E,expected_M\n";
    for (int i=0; i<burn_in_time.size(); i++){
        burn_in_csv << burn_in_time[i] << "," << expected_E_burn_in[i] 
                    << "," << expected_M_burn_in[i] << "\n";
    }
    burn_in_csv.close();
}


/**
 * @brief Produces sampled distribution of &epsilon; and m, and estimates for <&epsilon;>, <|m|>, C_v and &chi; 
 * 
 * @param N Numbers of MCMC-iterations
 * @param L Size of the Lattice (will have L * L elements)
 * @param T Temprature
 * @param burn_in_time Number of iterations to discard when producing estimates 
 */
void produce_distributions(const int iters, const int L, double T, int burn_in_time = 1000){ // note we don't have any indications on what the burn_in_time is yet
    IsingModel model(L, T, true);
    const int N = L * L;
    int sampled_energy[iters];
    int sampled_magnetization_abs[iters];
    for (int i = -burn_in_time; i < iters; i++){
        model.metropolis();
        if (i < 0) continue;
        sampled_energy[i] = model.get_energy();
        sampled_magnetization_abs[i] = abs(model.get_magnetization());
    }
    auto scale = [N](int x){return x / N;};
    map<int, double> buckets_epsilon = stat_utils::distribution(sampled_energy, iters, scale);
    map<int, double> buckets_m_abs = stat_utils::distribution(sampled_magnetization_abs, iters, scale);

    auto square = [](int x){return x * x;};
    double expected_epsilon = stat_utils::expected_value(buckets_epsilon);
    double expected_m_abs = stat_utils::expected_value(buckets_m_abs);
    double expected_energy = stat_utils::expected_value(sampled_energy, iters);
    double expected_energy_sq = stat_utils::expected_value(sampled_energy, iters, square);
    double expected_magnetization_abs = stat_utils::expected_value(sampled_magnetization_abs, iters);
    double expected_magnetization_sq = stat_utils::expected_value(sampled_magnetization_abs, iters, square);
    double c_v = (1. / N) * (1. / (T * T)) * (expected_energy_sq - expected_energy * expected_energy);
    double chi = (1. / N) * (1. / T) * (expected_magnetization_sq - expected_magnetization_abs * expected_magnetization_abs);

    cout << "Sampled distribution of epsilon" << endl;
    stat_utils::print_map(buckets_epsilon);
    cout << "Sampled distribution of |m|" << endl;
    stat_utils::print_map(buckets_m_abs);
    cout << "Other estimates: " << endl;
    cout << "<epsilon>: " <<  expected_epsilon << endl;
    cout << "<|m|>: " << expected_m_abs << endl;
    cout << "C_v: " <<  c_v << endl;
    cout << "Chi: " << chi << endl;
}



void test2x2(){
    const int N = 100000;
    int L = 2;
    double T = 1;
    produce_distributions(N, L, T);
}

int main(){
    //find_burn_in_time(200, 20, 1.0);
    // produce_distributions(10000, 10, 10.);
    test2x2();
}
