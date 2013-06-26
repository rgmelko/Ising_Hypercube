#ifndef PARAMSIM_H
#define PARAMSIM_H

#include <fstream>
#include <iostream>
using namespace std;

//Simple class to read in the simulation parameters from a file
class PARAMS
{
    public:
        int nX_;      //linear size of lattice
        int Dim_;      //dimension of lattice
        double Temp_;      //Temperature
        int EQL_;     //the number of Monte Carlo steps
        int MCS_;     //the number of Monte Carlo steps
        int nBin_;    //number of production bins
        long SEED_;

        PARAMS();
        void print();


}; //PARAMS

PARAMS::PARAMS(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");

    pfin >> nX_;
    pfin >> Dim_;
    pfin >> Temp_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> nBin_;
    pfin >> SEED_;
    pfin.close();

}//constructor


void PARAMS::print(){

    cout<<"Linear size "<<nX_<<endl;
    cout<<"Dimension "<<Dim_<<endl;
    cout<<"# Equil steps "<<EQL_<<endl;
    cout<<"# MC steps "<<MCS_<<endl;
    cout<<"# data bins "<<nBin_<<endl;
    cout<<"RNG seed "<<SEED_<<endl;

}

#endif
