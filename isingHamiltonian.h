#ifndef ISING HAMILTONIAN_H
#define ISING HAMILTONIAN_H

// isingHamiltonian.h
// a class calculate the energy of an Ising model

#include "spins.h"
#include "MersenneTwister.h"
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

        double Energy;  //total energy of the system

        //All the 2*D neighbors of a given site
        vector<vector<int> > All_Neighbors; 
        //you will double count if you calculate energy from this directly...

        IsingHamiltonian(Spins & sigma, HyperCube & cube); 
        double CalcEnergy(Spins & sigma);
        void LocalUpdate(Spins & sigma, double & T, MTRand & ran);

        void print();


};

//constructor
IsingHamiltonian::IsingHamiltonian(Spins & sigma, HyperCube & cube){

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


    cout<<CalcEnergy(sigma)/(1.0*N_)<<endl;      


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
double IsingHamiltonian::CalcEnergy(Spins & sigma){

    Energy = 0.0;

    for (int i=0; i<All_Neighbors.size(); i++){
        for (int j=0; j<All_Neighbors[i].size(); j++){
            Energy += -sigma.spin[i]*sigma.spin[All_Neighbors[i][j]];
        }//j
    }//i

    Energy /= 2.0; //double counting

    return Energy;

}

//Calculates a number of single-spin flips
void IsingHamiltonian::LocalUpdate(Spins & sigma, double & T, MTRand & ran){

    int site;  //random site for update
    double Ediff;
    double m_rand; //metropolis random number

    for (int j=0; j<N_; j++){ //peform N random single spin flips

        site = ran.randInt(N_-1);
        //cout<<"site is "<<site<<endl;

        Ediff = 0;
        for (int i=0; i<All_Neighbors[site].size(); i++)
            Ediff += -sigma.spin[site] * sigma.spin[All_Neighbors[site][i]];
        Ediff *= -2;

        //cout<<Energy<<" "<<Ediff<<endl;

        //Metropolis algorithm
        if (Ediff < 0){
            sigma.flip(site);
            Energy += Ediff;
        }
        else{
            m_rand = ran.rand();   // real number in [0,1]
            //cout<<"exponential "<<exp(-Ediff/T)<<" "<<m_rand<<endl;
            if ( exp(-Ediff/T) > m_rand){
                sigma.flip(site);
                Energy += Ediff;
            }
            // otherwise reject
            //else cout<<"reject: ";
        }

    }//j

    //cout<<"Emod "<<Energy<<endl;
}//LocalUpdate



#endif


