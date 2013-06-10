// A program simulate the Ising hamiltonian on a D-dimensional Hypercube
// Roger Melko, June 8, 2013
#include <iostream>
#include <vector>
using namespace std;

#include "hypercube.h"
#include "MersenneTwister.h"
#include "simparam.h"
#include "isingHamiltonian.h"
#include "measure.h"

int main(){

    double T; //temperature

    PARAMS param; //read parameter file
    //param.print();

    MTRand mrand(param.SEED_); //random number for metropolis

    HyperCube cube(param.nX_,param.Dim_); //initialize the lattice
    //cube.print();

    //define the Ising variables +1 or -1: initialize to 1
    Spins sigma(cube.N_);

    IsingHamiltonian hamil(sigma,cube);
    //hamil.print();

    //sigma.print();
    //sigma.flip(mrand.randInt(cube.N_-1));
    //cout<<"Energy: "<<hamil.CalcEnergy(sigma)<<endl;

    Measure accum(cube.N_,param);
    for (T = 20; T>0.4; T -= 0.4){
        for (int i=0; i<param.EQL_; i++) hamil.LocalUpdate(sigma,T,mrand);
        accum.zero();
        for (int i=0; i<param.MCS_; i++){ 
            hamil.LocalUpdate(sigma,T,mrand);
            accum.record(hamil.Energy,sigma);
        }
        accum.output(T);
    }

    return 0;

}
