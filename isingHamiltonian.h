#ifndef ISINGHAMILTONIAN_H
#define ISINGHAMILTONIAN_H

// isingHamiltonian.h
// a class calculate the energy of an Ising model
// note: BOOST_DISABLE_ASSERTS

#define PRINT_RED(x) std::cout << "\033[1;31m" << x << "\033[0m" << " "
#define PRINT_BLUE(x) std::cout << "\033[1;34m" << x << "\033[0m" << " "
#define PRINT_GREEN(x) std::cout << "\033[1;32m" << x << "\033[0m" << " "
#define PRINT_YELLOW(x) std::cout << "\033[1;33m" << x << "\033[0m" << " "

#include "spins.h"
#include "MersenneTwister.h"
#include <vector>
#include <iostream>
#include <boost/multi_array.hpp>

using namespace std;

//typedef vector<vector<int> > array_2t;
typedef boost::multi_array<int, 2> array_2t;
typedef boost::multi_array<int, 1> array_1t;

class IsingHamiltonian 
{
    public:
        int N_;   //number of lattice sites
        int D_;   //Dimension
        int L_;   //Linear size

        int Bonds_Per_Site; //total number of bonds per site

        double Energy;  //total energy of the system

        //All the 2*D neighbors of a given site
        array_2t All_Neighbors; 
		//the above also applies for usual Fortuin Kasteleyn clusters
		array_1t occupancy;

		//The neighbors from a percolation perspective

        IsingHamiltonian(Spins & sigma, HyperCube & cube); 
		void print();
        double CalcEnergy(Spins & sigma);
        void LocalUpdate(Spins & sigma, double & T, MTRand & ran);
        void CalculateOccupancy(Spins & sigma);

};

//constructor
IsingHamiltonian::IsingHamiltonian(Spins & sigma, HyperCube & cube){

    L_ = cube.L_;
    D_ = cube.D_;
    N_ = cube.N_;

    sigma.resize(N_); //these are the degrees of freedom (1 cells)
    sigma.randomize();

	occupancy.resize(boost::extents[N_]); //calculate simple site percolation
	for (int j=0; j<sigma.spin.size(); j++)
		occupancy[j] = (sigma.spin[j]+1)/2; //0 or 1

    Bonds_Per_Site = 2*D_;  //this will double count the total number of bonds  

    //resize the empty 2D array
    All_Neighbors.resize(boost::extents[N_][Bonds_Per_Site]);

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

    //PRINT_GREEN("L, D and N: ");
    //cout<<L_<<" "<<D_<<" "<<N_<<endl;

    //PRINT_BLUE("Boost array dimension and shape: ");
    //cout<<All_Neighbors.dimensionality<<" ";
    //cout<<All_Neighbors.shape()[0]<<" ";
    //cout<<All_Neighbors.shape()[1]<<endl;

    //for (int i=0; i<All_Neighbors.size(); i++){
    //    cout<<i<<" ";
    //    for (int j=0; j<Bonds_Per_Site; j++)
    //        cout<<All_Neighbors[i][j]<<" ";
    //    cout<<endl;
    //}//i

    PRINT_RED("occupancy");cout<<endl;
	for (int i=0; i<N_; i++){
		cout<<occupancy[i]<<" ";
		//if ((i+1)%L_ == 0) cout<<endl;
	}
	cout<<endl;

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


//This calculates the occupancy, of arbitrary definition, of clusters for the 
// percolation calculation.
void IsingHamiltonian:: CalculateOccupancy(Spins & sigma){

    for (int j=0; j<sigma.spin.size(); j++)
        occupancy[j] = (sigma.spin[j]+1)/2; //0 or 1

}


#endif


