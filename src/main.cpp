#include "project4/ising_model.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main(){
    int L = 2;
    IsingModel model(L);
    vector<vector<int>> spins = model.get_spins();
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            cout << spins[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "energy: " << model.get_energy() <<  "\n";
    model.metropolis();
}
