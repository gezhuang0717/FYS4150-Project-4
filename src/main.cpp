#include "project4/ising_model.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

map<int, float> probability_distrubution(int *samples, int N){
    map<int, float> buckets;
    for (int i = 0; i < N; i++){
        int sample = samples[i];
        buckets.emplace(sample, 0);
        buckets[sample] += 1./N;
    }
    return buckets;
}

int test2x2(){
    const int N = 100000;
    int L = 2;
    double T = 10;
    IsingModel model(L, T);
    int sampled_E[N];
    int sampled_M[N];
    for (int i = -10000; i < N; i++){
        model.metropolis();
        if (i < 0) continue;
        sampled_E[i] = model.get_energy();
        sampled_M[i] = model.get_magnetisation();
    }
    map<int, float> buckets_E = probability_distrubution(sampled_E, N);
    map<int, float> buckets_M = probability_distrubution(sampled_M, N);
    for (auto const& bucket : buckets_E){
        cout << bucket.first << ": " << bucket.second << endl;
    }
    cout << "------" << endl;
    float E_M = 0;
    for (auto const& bucket : buckets_M){
        cout << bucket.first << ": " << bucket.second << endl;
        E_M += bucket.first * bucket.second;
    }
    cout << E_M << endl;
}

int main(){
    test2x2();
}
