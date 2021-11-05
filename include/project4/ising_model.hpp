#pragma once

#include <vector>
#include <random>
#include <iostream>

using namespace std;
class IsingModel{
    public:
        IsingModel(int L);
        vector<vector<int>> get_spins();
        // Performs one Monte Carlo cycle of Metropolis algorithm and
        // Updates energy- and magnetization-values.
        void metropolis();
        // Should only be called by constructor, but kept public as it might be useful
        // for debugging later.
        void set_energy();
        int get_energy();
    private:
        // Called by constructor to initialise 2D vector of spins.
        // (currently a placeholder-function that sets every other spin equal to -1)
        void initialise_spins(int L);
        mt19937 rng;
        uniform_int_distribution<int> rand_index;
        vector<vector<int>> spins;
        int L;
        int E;
};
