#include "project4/ising_model.hpp"
#include "project4/map_utils.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;


void produce_distributions(const int N, int L, double T, int burn_in_time = 1000){ // note we don't have any indications on what the burn_in_time is yet
    IsingModel model(L, T);
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
    produce_distributions(10000, 10, 10.);
}
