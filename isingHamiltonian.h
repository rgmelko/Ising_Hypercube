#ifndef ISING HAMILTONIAN_H
#define ISING HAMILTONIAN_H

// isingHamiltonian.h
// a class calculate the energy of an Ising model

#include <vector>
#include <iostream>

using namespace std;

class IsingHamiltonian 
{
    public:
        int N_;   //number of lattice sites
        int D_;   //Dimension
        int L_;   //Linear size

        int Bonds_Per_Site; //total number of bonds per site

        int Energy;  //total energy of the system

        //All the 2*D neighbors of a given site
        vector<vector<int> > All_Neighbors; 
        //you will double count if you calculate energy from this directly...

        IsingHamiltonian(vector <int> &, HyperCube & cube); 
        //double CalcEnergy(vector <int> &);

        void print();


};

//constructor
IsingHamiltonian::IsingHamiltonian(vector <int> & Spins, HyperCube & cube){

    L_ = cube.L_;
    D_ = cube.D_;
    N_ = cube.N_;

    Bonds_Per_Site = 2*D_;  //this will double count the total number of bonds  

    //resize the empty 2D array
    All_Neighbors.resize(N_);
    for (int i=0; i<All_Neighbors.size(); i++)
        All_Neighbors[i].resize(Bonds_Per_Site);

    //build it from the hypercubic lattice
    for (int i=0; i<All_Neighbors.size(); i++){
        for (int j=0; j<Bonds_Per_Site/2; j++){
            All_Neighbors[i][j]= cube.Neighbors[i][j];
            All_Neighbors[cube.Neighbors[i][j]][j+Bonds_Per_Site/2]= i;
        }//j
    }//i


    //CalcEnergy(Spins);      


}//constructor


//print
void IsingHamiltonian::print(){

    cout<<L_<<" "<<D_<<" "<<N_<<endl;

    for (int i=0; i<All_Neighbors.size(); i++){
        cout<<i<<" ";
        for (int j=0; j<Bonds_Per_Site; j++)
            cout<<All_Neighbors[i][j]<<" ";
        cout<<endl;
    }//i

}//print


//loops through to calculate the energy
//IsingHamiltonian::CalcEnergy(vector<int> & Spins){
//
//    for (int i=0; i< N_; i++){
//
//
//    }//i
//
//}


#endif