#include "project4/ising_model.hpp"

IsingModel::IsingModel(int lattice_length, double T){
    L = lattice_length;
    beta = 1/T; 
    int seed = 7773;
    rng = mt19937(seed);
    rand_index = uniform_int_distribution<int>(0, L-1);
    uniform = uniform_real_distribution<double>(0, 1);
    initialize_spins(L);
    set_energy();
    set_magnetisation();
}

// Placeholder initialisations function.
void IsingModel::initialize_spins(int L){
    vector<vector<int>> initial(L, vector<int>(L, 1));
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            if ((i+j)%2 == 0){
                initial[i][j] = -1 * initial[i][j];
            }
        }
    }
    spins = initial;
}

void IsingModel::set_energy(){
    int energy = 0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            energy += spins[i][j] * spins[(i+1)%L][j];
            energy += spins[i][j] * spins[i][(j+1)%L];
        }
    }
    E = energy;
}

void IsingModel::set_magnetisation(){
    int magnetisation = 0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            magnetisation += spins[i][j];
        }
    }
    M = magnetisation;
}

int IsingModel::get_magnetisation(){
    return M;
}

int IsingModel::get_energy(){
    return E;
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
        delta_E = -2 * spins[ix][iy] * (spins[ix][(iy + 1) % L]
                                    + spins[ix][(iy - 1 + L) % L]
                                    + spins[(ix + 1) % L][iy]
                                    + spins[(ix - 1 + L) % L][iy]);
        //compute exp(-beta*delta_E)
        // TODO: replace with pre-computed lookup
        double w = exp(-beta * (double)delta_E);
        double r = uniform(rng);
        // Check if spin should be flipped
        if (delta_E < 0 or r <= w){
            // Flip spin
            spins[ix][iy] *= -1;
            // Update energy
            E = delta_E;
            // Update magnetisation
            M += 2 * spins[ix][iy];
        }
    }

}
