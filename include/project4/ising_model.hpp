#pragma once

#include <vector>
#include <random>
#include <iostream>

using namespace std;
class IsingModel{
    public:
        IsingModel(int L);
        vector<vector<int>> get_spins();
        void metropolis();
        void set_energy();
        int get_energy();
    private:
        void initialise_spins(int L);
        mt19937 rng;
        uniform_int_distribution<int> rand_index;
        vector<vector<int>> spins;
        int L;
        int E;
};