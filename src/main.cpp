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
    const int N = 100;
    int L = 2;
    double T = 10;
    IsingModel model(L, T);
    int sampled_E[N];
    for (int i = 0; i < N; i++){
        model.metropolis();
        sampled_E[i] = model.get_energy();
    }
    map<int, float> buckets = probability_distrubution(sampled_E, N);
    for (auto const& bucket : buckets){
        cout << bucket.first << ": " << bucket.second << endl;
    }
}

int main(){
    int L = 2;
    double T = 10;
    IsingModel model(L, T);
    vector<vector<int>> spins = model.get_spins();
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            cout << spins[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "energy: " << model.get_energy() <<  "\n";

    model.metropolis();
    model.metropolis();
    model.metropolis();

    spins = model.get_spins();
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            cout << spins[i][j] << " ";
        }
        cout << "\n";
    }

    test2x2();
}
