// A program to simulate the Ising hamiltonian on a D-dimensional Hypercube
// Roger Melko, June 8, 2013
#include <iostream>
#include <vector>
using namespace std;

#include "hypercube.h"
#include "MersenneTwister.h"
#include "simparam.h"
//#include "generalD_1_2.code.h"
#include "isingHamiltonian.h"
#include "measure.h"

typedef boost::multi_array<int, 2> array_2t;

int main(){

    PARAMS param; //read parameter file
    //param.print();

    MTRand mrand(param.SEED_); //random number for metropolis

    HyperCube cube(param.nX_,param.Dim_); //initialize the lattice
    //cube.print();

    //define the Ising variables +1 or -1: initialize to 1
    //Spins sigma(cube.N_);
    Spins sigma; //Assign size of spins in Hamiltonian below

    IsingHamiltonian hamil(sigma,cube);
    //GeneralD12Code hamil(sigma,cube);
    hamil.print();

    return 0;

    Measure accum(hamil.N_,param);     //Ising model
    //Measure accum(hamil.N1,param);  //toric code

    //insert T loop here
    for (double T = param.Temp_; T>0.5; T-=0.1){

        //Equilibriation
        for (int i=0; i<param.EQL_; i++) hamil.LocalUpdate(sigma,T,mrand);

        //MCS binning
        for (int k=0; k<param.nBin_; k++){ 
            accum.zero();
            for (int i=0; i<param.MCS_; i++){ 
                hamil.LocalUpdate(sigma,T,mrand);
                accum.record(hamil.Energy,sigma);
            }//i
            accum.output(T);
            //sigma.print();
        }//k
    }//T

    return 0;

}
