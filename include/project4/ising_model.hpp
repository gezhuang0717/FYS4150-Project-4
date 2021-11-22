#pragma once

#include <map>
#include <vector>
#include <random>
#include <iostream>

using namespace std;
class IsingModel{
    public:
        IsingModel(int L, double T, bool random_spins, int seed);
        /**
         * @brief Get the matrix of spins
         * 
         * @return vector<vector<int>>
         */
        vector<vector<int>> get_spins();
        /**
         * @brief Performs one iteration of Metropolis algorithm
         * 
         */
        void metropolis();
        /**
         * @brief Get the energy of the IsingModel
         * 
         * @return int 
         */
        int get_energy();
        /**
         * @brief Get the magnetization of the IsingModel
         * 
         * @return int 
         */
        int get_magnetization();
        /**
         * @brief Get m of the IsingModel
         * 
         * @return double 
         */
        double get_m();
        /**
         * @brief Get w from the precomputed set of exponential values
         * 
         * @param delta_E The change of energy with this spin flip
         * @return double 
         */
        double get_w(int delta_E);
        /**
         * @brief Get the epsilon of the IsingModel
         * 
         * @return double 
         */
        double get_epsilon();
        /**
         * @brief Print the spin state
         * 
         */
        void print();
    private:
        
        /**
         * @brief 
         * Called by constructor to initialise 2D vector of spins.
         * 
         * @param L Size of the ising model
         */
        void initialize_spins(int L);
        /**
         * @brief Set the initial magnetization
         * 
         */
        void set_magnetization();
        /**
         * @brief Set the initial energy
         * 
         */
        void set_energy();
        mt19937 rng;
        uniform_int_distribution<int> rand_index;
        uniform_int_distribution<int> rand_1_or_0;
        uniform_real_distribution<double> uniform;
        vector<vector<int>> spins;
        int L;
        int E;
        int M;
        double beta;
        bool rand_spins;
        /**
         * @brief Precomputes exponential values that are needed every MC-cycle, but are limited to five possible values
         * 
         */
        void precompute_exp_factors();
        double exp_factors[5];
};
