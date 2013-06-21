#ifndef SPINS_H
#define SPINS_H

// spins.h
// a small class that contains a vector of lattice Ising spins

#include <vector>
#include <iostream>

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

//a single-spin flip
void Spins::flip(int index){

    spin.at(index) *= -1;

}//flip

//a print function
void Spins::print(){

    for (int i=0;i<spin.size();i++){
        cout<<spin[i]<<" ";
    }//i
    cout<<endl;

}//print

#endif
