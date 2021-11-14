#include "project4/ising_model.hpp"

IsingModel::IsingModel(int lattice_length, double T, bool random_spins){
    L = lattice_length;
    beta = 1/T; 
    rand_spins = random_spins;
    // int seed = 7773;
    // rng = mt19937(seed);
    rng = mt19937(69);
    rand_index = uniform_int_distribution<int>(0, L-1);
    rand_1_or_0 = uniform_int_distribution<int>(0, 1);
    uniform = uniform_real_distribution<double>(0, 1);
    initialize_spins(L);
    set_energy();
    set_magnetization();
}

// Placeholder initialisations function.
void IsingModel::initialize_spins(int L){
    vector<vector<int>> initial(L, vector<int>(L, 1));
    if (rand_spins){
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                initial[i][j] = initial[i][j] - 2*rand_1_or_0(rng);
            }
        }
    }
    spins = initial;
}

void IsingModel::set_energy(){
    int energy = 0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            energy += - spins[i][j] * spins[(i+1)%L][j];
            energy += - spins[i][j] * spins[i][(j+1)%L];
        }
    }
    E = energy;
}

void IsingModel::set_magnetization(){
    int magnetization = 0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            magnetization += spins[i][j];
        }
    }
    M = magnetization;
}

int IsingModel::get_magnetization(){
    return M;
}

int IsingModel::get_energy(){
    return E;
}

float IsingModel::get_epsilon(){
    return (float)E / (L * L);
}

float IsingModel::get_m(){
    return (float)M / (L * L);
}

vector<vector<int>> IsingModel::get_spins(){
    return spins;
}


void IsingModel::metropolis(){
    int ix, iy, delta_E;
    for (int _ = 0; _ < L * L; _++){
        // Draw 2 random indices
        ix = rand_index(rng);
        iy = rand_index(rng);
        // Compute energy difference for spin flip
        delta_E = 2 * spins[ix][iy] * (spins[ix][(iy + 1) % L]
                                    + spins[ix][(iy - 1 + L) % L]
                                    + spins[(ix + 1) % L][iy]
                                    + spins[(ix - 1 + L) % L][iy]);
        //compute exp(-beta*delta_E)
        // TODO: replace with pre-computed lookup
        double w = exp(-beta * (double)delta_E);
        double r = uniform(rng);
        // Check if spin should be flipped
        if (r <= w){
            // Flip spin
            spins[ix][iy] *= -1;
            // Update energy
            E += delta_E;
            // Update magnetization
            M += 2 * spins[ix][iy];
        }
    }
}

void IsingModel::print(){
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            cout << spins[i][j] << " ";
        }
        cout << "\n";
    }
}