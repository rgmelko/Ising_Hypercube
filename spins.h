#ifndef SPINS_H
#define SPINS_H

// spins.h
// a small class that contains a vector of lattice Ising spins

#include <vector>
#include <iostream>
#include "MersenneTwister.h"

using namespace std;

class Spins
{
    public:
        int N_; //total number of lattice sites

        //the lattice is a vector of vectors: no double counting
        vector<int> spin;

        //public functions
        Spins(int N);
        Spins();
        void resize(int N);
        void flip(int index);
        void print();
        void randomize();

};

//constructor 1
//takes the total number of lattice sites
Spins::Spins(){

    spin.clear(); 

}

//constructor 2
//takes the total number of lattice sites
Spins::Spins(int N){

    N_ = N;

    spin.resize(N_,1); //assign every spin as 1

}

//takes the total number of lattice sites
void Spins::resize(int N){

    N_ = N;

    spin.resize(N_,1); //assign every spin as 1

}

void Spins::randomize(){

    MTRand irand(129345); //random number 

    int ising_spin;
    for (int i = 0; i<spin.size(); i++){
        ising_spin = 2*irand.randInt(1)-1;
        //cout<<ising_spin<<" ";
        spin.at(i) = ising_spin;
    }


}//randomize

//a single-spin flip
void Spins::flip(int index){

    spin.at(index) *= -1;

}//flip


//a print function
void Spins::print(){

    for (int i=0;i<spin.size();i++){
        cout<<(spin[i]+1)/2<<" ";
    }//i
    cout<<endl;

}//print

#endif
