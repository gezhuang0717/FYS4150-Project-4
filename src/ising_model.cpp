#include "project4/ising_model.hpp"

IsingModel::IsingModel(int lattice_length){
    int seed = 7773;
    rng = mt19937(seed);
    L = lattice_length;
    rand_index = uniform_int_distribution<int>(0, L-1);
    initialise_spins(L);
    set_energy();
}

// Placeholder initialise function.
void IsingModel::initialise_spins(int L){
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

int IsingModel::get_energy(){
    return E;
}

vector<vector<int>> IsingModel::get_spins(){
    return spins;
}

void IsingModel::metropolis(){
    int ix, iy, delta_E;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            // Draw 2 random indices
            ix = rand_index(rng);
            iy = rand_index(rng);

            // Compute energy difference for spin flip

            // Får segmentation fault her pga. % er remainder-operator
            // og ikke modulo, så får negative indekser. Skal fikse senere.
            delta_E = 2*spins[ix][iy] *(spins[ix][(iy+1)%L]
                                        + spins[ix][(iy-1)%L]
                                        + spins[(ix+1)%L][iy]
                                        + spins[(ix-1)%L][iy]);

            // To be implemented:

            // Check if spin should be flipped
            if (...){
                // Flip spin
                ...
                // Update energy
            }
        }
        
    }

}