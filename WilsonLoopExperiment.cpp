// A program to simulate the Ising and (1,d-1) codes on a D-dimensional Hypercube
// Roger Melko, June 8, 2013
//
// Requires BOOST multi_array: http://www.boost.org
// compile example:  g++ -O3 IsingMonteCarlo.cpp -I /opt/local/include/  
//
#include <iostream>
#include <vector>
using namespace std;

#include <boost/multi_array.hpp>

#include "hypercube.h"
#include "MersenneTwister.h"
#include "simparam.h"
#include "generalD_1_2.code.h"
//#include "isingHamiltonian.h"
#include "measure.h"
#include "percolation.h"

int main ( int argc, char *argv[] )
{
    //A different command line argument for each simulation
    int seed_add = 0;
    //if ( argc != 2 ){ 
    //    cout<<"usage: "<< argv[0] <<" integer \n";
    //    return 1;
    //}
    //else {
    //    seed_add = strtol(argv[1], NULL, 10);
    //}

    //First, we call several constructors for the various objects used

    PARAMS param; //read parameter file: L, D, T, etc.  See param.data

    MTRand mrand(param.SEED_+seed_add); //random number generator

    HyperCube cube(param.nX_,param.Dim_); //initialize the lattice

    //define the Ising variables +1 or -1 
    Spins sigma; //Assign number of spins in the Hamiltonian below

    GeneralD12Code hamil(sigma,cube,param.H_); //toric code
	//--
	//sigma.flip(1);
	//cout<<hamil.CalcEnergyDiff(sigma,1,param.H_)<<" ";
	//cout<<hamil.CalcEnergy(sigma,param.H_)<<endl;
	//sigma.flip(6);
	//cout<<hamil.CalcEnergyDiff(sigma,1,param.H_)<<" ";
	//cout<<hamil.CalcEnergy(sigma,param.H_)<<endl;

    //hamil.PreparePercolation(sigma,cube); //for D>2 toric code percolation only

    Measure accum(hamil.N1,param);  //toric code

    double H = param.H_;
    double T = param.Temp_;
    //for (T = param.Temp_; T<param.Tlow_; T+=param.Tstep_){ //up
        //Equilibriation
        for (int i=0; i<param.EQL_; i++) {
            hamil.LocalUpdate(sigma,T,mrand,H);
        }
        //MCS binning
        for (int k=0; k<param.nBin_; k++){ 
            accum.zero();
            for (int i=0; i<param.MCS_; i++){ 
                hamil.LocalUpdate(sigma,T,mrand,H);
                accum.record(hamil.Energy,sigma,hamil.WilsonLoops);
                //accum.outputWilsonLoop(sigma,hamil.WilsonLoops,seed_add);

            }//i
            accum.output(T,H);
        }//k
    //}//T up

//    //Stationary...
//    //Equilibriation
//    //for (int i=0; i<param.EQL_; i++) {
//    //    hamil.LocalUpdate(sigma,T,mrand);
//    //}
//    //MCS binning
//    for (int k=0; k<param.nBin_; k++){ 
//        accum.zero();
//        for (int i=0; i<param.MCS_; i++){ 
//            hamil.LocalUpdate(sigma,T,mrand);
//            accum.record(hamil.Energy,sigma,hamil.WilsonLoops);
//            accum.outputWilsonLoop(sigma,hamil.WilsonLoops,seed_add);
//
//        }//i
//        accum.output(T);
//    }//k
//
//    //down
//    for (T = param.Tlow_; T>param.Temp_; T-=param.Tstep_){ //up
//        //Equilibriation
//        //for (int i=0; i<param.EQL_; i++) {
//        //    hamil.LocalUpdate(sigma,T,mrand);
//        //}
//        //MCS binning
//        for (int k=0; k<param.nBin_; k++){ 
//            accum.zero();
//            for (int i=0; i<param.MCS_; i++){ 
//                hamil.LocalUpdate(sigma,T,mrand);
//                accum.record(hamil.Energy,sigma,hamil.WilsonLoops);
//                //accum.outputWilsonLoop(sigma,hamil.WilsonLoops,seed_add);
//
//            }//i
//            accum.output(T);
//        }//k
//    }//T up

    //accum.outputWilsonLoop(sigma,hamil.WilsonLoops,seed_add);

    return 0;

}
