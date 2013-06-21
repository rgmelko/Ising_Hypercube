// A program simulate the Ising hamiltonian on a D-dimensional Hypercube
// Roger Melko, June 8, 2013
#include <iostream>
#include <vector>
using namespace std;

#include "hypercube.h"
#include "MersenneTwister.h"
#include "simparam.h"
#include "threeD_1_2_code.h"
#include "measure.h"

int main(){

    double T; //temperature

    PARAMS param; //read parameter file
    //param.print();

    MTRand mrand(param.SEED_); //random number for metropolis

    HyperCube cube(param.nX_,param.Dim_); //initialize the lattice
    //cube.print();

    //define the Ising variables +1 or -1: initialize to 1
    //Spins sigma(cube.N_);
	Spins sigma; //Assign size of spins in Hamiltonian

	ThreeD12Code hamil(sigma,cube);
	hamil.print();

	return 1;

    //sigma.print();
    //sigma.flip(mrand.randInt(cube.N_-1));
    //cout<<"Energy: "<<hamil.CalcEnergy(sigma)<<endl;
    //sigma.print();

    Measure accum(cube.N_,param);
    for (T = 4; T>0.1; T -= 0.1){
        for (int i=0; i<param.EQL_; i++) hamil.LocalUpdate(sigma,T,mrand);
        accum.zero();
        for (int i=0; i<param.MCS_; i++){ 
            hamil.LocalUpdate(sigma,T,mrand);
            accum.record(hamil.Energy,sigma);
        }
        accum.output(T);
    }

    sigma.print();

    return 0;

}
