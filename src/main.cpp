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
        map<int, float> buckets_E = stat_utils::distribution(sampled_E + burn_in, N - burn_in);
        map<int, float> buckets_M = stat_utils::distribution(sampled_M + burn_in, N - burn_in);
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
 * @brief Produces sampled distribution of &epsilon; and m, and estimates for <&eplsilon;>, <|m|>, C_v and &chi 
 * 
 * @param N Numbers of MCMC-iterations
 * @param L Size of the Lattice (will have L * L elements)
 * @param T Temprature
 * @param burn_in_time Number of iterations to discard when producing estimates 
 */
void produce_distributions(const int N, const int L, double T, int burn_in_time = 1000){ // note we don't have any indications on what the burn_in_time is yet
    IsingModel model(L, T, true);
    int sampled_energy[N];
    int sampled_magnetization_abs[N];
    for (int i = -burn_in_time; i < N; i++){
        model.metropolis();
        if (i < 0) continue;
        sampled_energy[i] = model.get_energy();
        sampled_magnetization_abs[i] = abs(model.get_magnetization());
    }
    auto scale = [L](int x){return x / (L * L);};
    map<int, float> buckets_epsilon = stat_utils::distribution(sampled_energy, N, scale);
    map<int, float> buckets_m_abs = stat_utils::distribution(sampled_magnetization_abs, N, scale);

    auto square = [L](int x){return x * x;};
    float expected_epsilon = stat_utils::expected_value(buckets_epsilon);
    float expected_m_abs = stat_utils::expected_value(buckets_m_abs);
    float expected_energy = stat_utils::expected_value(sampled_energy, N);
    float expected_energy_sq = stat_utils::expected_value(sampled_energy, N, square);
    float expected_magnetization_abs = stat_utils::expected_value(sampled_magnetization_abs, N);
    float expected_magnetization_sq = stat_utils::expected_value(sampled_magnetization_abs, N, square);

    cout << "Sampled distribution of epsilon" << endl;
    stat_utils::print_map(buckets_epsilon);
    cout << "Sampled distribution of |m|" << endl;
    stat_utils::print_map(buckets_m_abs);
    cout << "Other estimates: " << endl;
    cout << "<epsilon>: " <<  expected_epsilon << endl;
    cout << "<|m|>: " << expected_m_abs << endl;
    cout << "C_v: " << (1. / N) * (1. / (T * T)) * (expected_energy_sq - expected_energy * expected_energy) << endl;
    cout << "Chi: " << (1. / N) * (1. / T) * (expected_magnetization_sq - expected_magnetization_abs * expected_magnetization_abs) << endl;
}



void test2x2(){
    const int N = 100000;
    int L = 2;
    double T = 10;
    produce_distributions(N, L, T);
}

int main(){
    //find_burn_in_time(200, 20, 1.0);
    // produce_distributions(10000, 10, 10.);
    test2x2();
}
