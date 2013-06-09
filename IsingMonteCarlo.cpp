// A program simulate the Ising hamiltonian on a D-dimensional Hypercube
// Roger Melko, June 8, 2013
#include <iostream>
#include <vector>
using namespace std;

#include "hypercube.h"
#include "MersenneTwister.h"
#include "simparam.h"
#include "isingHamiltonian.h"

int main(){

    double T=1; //temperature

    PARAMS param; //read parameter file
    //param.print();

    MTRand mrand(param.SEED_); //random number for metropolis

    HyperCube cube(param.nX_,param.Dim_); //initialize the lattice
    //cube.print();

    //define the Ising variables +1 or -1: initialize to 1
    Spins sigma(cube.N_);

    IsingHamiltonian hamil(sigma,cube,T);
    //hamil.print();

	sigma.flip(0);
	cout<<hamil.CalcEnergy(sigma,T)<<endl;

    return 0;

}
