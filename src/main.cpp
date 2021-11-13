#include "project4/ising_model.hpp"
#include "project4/map_utils.hpp"

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
        sampled_M[i] = model.get_magnetisation();
        model.metropolis();
    }
    vector<int> burn_in_time;
    vector<double> expected_E_burn_in;
    vector<double> expected_M_burn_in;
    for (int burn_in=0; burn_in<N; burn_in+=max(1, N/100)){
        map<int, float> buckets_E = map_utils::distribution(sampled_E + burn_in, N - burn_in);
        map<int, float> buckets_M = map_utils::distribution(sampled_M + burn_in, N - burn_in);
        expected_E_burn_in.push_back(map_utils::expected_value(buckets_E));
        expected_M_burn_in.push_back(map_utils::expected_value(buckets_M));
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

void produce_distributions(const int N, int L, double T, int burn_in_time = 1000){ // note we don't have any indications on what the burn_in_time is yet
    IsingModel model(L, T, true);
    int sampled_E[N];
    int sampled_M[N];
    for (int i = -burn_in_time; i < N; i++){
        model.metropolis();
        if (i < 0) continue;
        sampled_E[i] = model.get_energy();
        sampled_M[i] = model.get_magnetisation();
    }
    map<int, float> buckets_E = map_utils::distribution(sampled_E, N);
    map<int, float> buckets_M = map_utils::distribution(sampled_M, N);
    cout << "Sampled distribution of E" << endl;
    map_utils::print_map(buckets_E);
    cout << "Sampled expected E: " << map_utils::expected_value(buckets_E) << endl;
    cout << "------" << endl;
    cout << "Sampled distribution of M" << endl;
    float E_M = map_utils::expected_value(buckets_M);
    map_utils::print_map(buckets_M);
    cout << "Sampled expected M: " << E_M << endl;
}



void test2x2(){
    const int N = 100000;
    int L = 2;
    double T = 10;
    produce_distributions(N, L, T);
}

int main(){
    find_burn_in_time(200, 20, 1.0);
    // produce_distributions(10000, 10, 10.);
}
