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
    for (int i = -10000; i < N; i++){
        model.metropolis();
        if (i < 0) continue;
        sampled_E[i] = model.get_energy();
    }
    map<int, float> buckets = probability_distrubution(sampled_E, N);
    for (auto const& bucket : buckets){
        cout << bucket.first << ": " << bucket.second << endl;
    }
}

int main(){
    test2x2();
}
