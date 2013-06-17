#ifndef TWODTORICCODE_H
#define TWODTORICCODE_H

// 2DtoricCode.h
// a class to perform a simple metropolis MC on a 2D Sz toric code

#include "spins.h"
#include "MersenneTwister.h"
#include <vector>
#include <iostream>

using namespace std;

class TwoDToricCode
{
    public:
        int N_;   //number of lattice sites
        int D_;   //Dimension
        int L_;   //Linear size

        double Energy;  //total energy of the system

        //All the 2*D neighbors of a given site
        vector<vector<int> > All_Neighbors; 

        //The sigma-z plaquette
        vector<vector<int> > Plaquette;

        TwoDToricCode(Spins & sigma, HyperCube & cube); 
        double CalcEnergy(Spins & sigma);
        void LocalUpdate(Spins & sigma, double & T, MTRand & ran);

        void print();

    private:

        int Bonds_Per_Site; //total number of bonds per site

};

//constructor
TwoDToricCode::TwoDToricCode(Spins & sigma, HyperCube & cube){

    L_ = cube.L_;
    D_ = 2;
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

    //use it to built the sigma-z plaquettes
    vector <int> temp;
    temp.assign(4,0);  //assign 4 zeros to this vector
    for (int i=0; i<N_; i++){
        temp[0] = i;
        temp[1] = All_Neighbors[i][0];
        temp[2] = All_Neighbors[i][1];
        temp[3] = All_Neighbors[temp[1]][1];
        Plaquette.push_back(temp);
    }//i


    cout<<CalcEnergy(sigma)<<endl;      
    //cout<<CalcEnergy(sigma)/(1.0*N_)<<endl;      


}//constructor


//print
void TwoDToricCode::print(){

    cout<<L_<<" "<<D_<<" "<<N_<<endl;

    for (int i=0; i<All_Neighbors.size(); i++){
        cout<<i<<" ";
        for (int j=0; j<Bonds_Per_Site; j++)
            cout<<All_Neighbors[i][j]<<" ";
        cout<<endl;
    }//i

    cout<<"Plaquette \n";
    for (int i=0; i<Plaquette.size(); i++){
        cout<<i<<" ";
        for (int j=0; j<4; j++)
            cout<<Plaquette[i][j]<<" ";
        cout<<endl;
    }//i

}//print


//loops through to calculate the energy
double TwoDToricCode::CalcEnergy(Spins & sigma){

    double eTemp = 0.0;

    for (int i=0; i<Plaquette.size(); i++){
        eTemp -=  sigma.spin[Plaquette[i][0]]*sigma.spin[Plaquette[i][1]]
                  *sigma.spin[Plaquette[i][2]]*sigma.spin[Plaquette[i][3]];
    }//i

    return eTemp;

}

//Calculates a number of single-spin flips
void TwoDToricCode::LocalUpdate(Spins & sigma, double & T, MTRand & ran){

    int site;  //random site for update
    double Eold, Enew, Ediff;
    double m_rand; //metropolis random number

    for (int j=0; j<N_; j++){ //peform N random single spin flips

        site = ran.randInt(N_-1);
        //cout<<"site is "<<site<<endl;

        sigma.flip(site);  //trial flip
        Eold = Energy;
        Enew = CalcEnergy(sigma);
        Ediff = Enew - Eold;

        //cout<<Energy<<" "<<Ediff<<endl;

        //Metropolis algorithm
        if (Ediff < 0){
            Energy = Enew;
        }
        else{
            m_rand = ran.rand();   // real number in [0,1]
            //cout<<"exponential "<<exp(-Ediff/T)<<" "<<m_rand<<endl;
            if ( exp(-Ediff/T) > m_rand){
                Energy = Enew;
            }
            else{ // otherwise reject
                sigma.flip(site);
                Energy = Eold; //redundant
            }
        }

    }//j

    //cout<<"Emod "<<Energy<<endl;
}//LocalUpdate



#endif


