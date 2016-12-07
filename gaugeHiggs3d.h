#ifndef THREED_1_2_CODE_H
#define THREED_1_2_CODE_H

// threeD_1_2_code.h
// a class to perform a simple metropolis MC on a 2D (1,2) toric code

#include "spins.h"
#include "MersenneTwister.h"
#include <vector>
#include <iostream>

using namespace std;

class ThreeD12Code
{
    public:
        int N1;   //number of DEGREES OF FREEDOM (1-cells)
        int D_;   //Dimension
        int L_;   //Linear size

        double Energy;  //total energy of the system

        //All the 1-cells (bonds) that are attached to 2-cells (faces)
        vector<vector<int> > All_Neighbors; 

        //The Face operators
        vector<vector<int> > Plaquette;
        ThreeD12Code(Spins & sigma, HyperCube & cube); 
        double CalcEnergy(Spins & sigma);
        double CalcEnergyDiff(Spins & sigma, const int & flipsite);

        void LocalUpdate(Spins & sigma, const double & T, MTRand & ran);

        void print();

    private:

        int Faces; //The total number of faces

};

//constructor
ThreeD12Code::ThreeD12Code(Spins & sigma, HyperCube & cube){

    L_ = cube.L_;
    D_ = 3;
    N1 = 3*cube.N_; //3D lattice BOND variables...

    sigma.resize(N1); //these are the degrees of freedom
    sigma.randomize();

    Faces = N1;  //in 3D the number of faces equals the numbers of DOFs

    //use it to built the sigma-z plaquettes
    vector <int> temp;
    temp.assign(4,0);  //assign 4 zeros to this vector

    int Xneigh;
    int Yneigh;
    int Zneigh;
    for (int i=0; i<Faces; i += 3 ){
        //XY plane
        temp[0] = i;
        temp[1] = i+1;
        temp[2] = i+4;
        temp[3] = i+3*L_;
        //fix boundaries
        if ((i/3+1)%L_ == 0)
            temp[2] -= 3*L_;
        if ( ( (i/3)%(L_*L_) >= (L_*L_) - L_) && ( (i/3)%(L_*L_) < (L_*L_) ) )
            temp[3] -= 3*L_*L_;

        Plaquette.push_back(temp);

        //XZ plane
        temp[0] = i;
        temp[1] = i+2;
        temp[2] = i+5;
        temp[3] = i+3*L_*L_;
        if ((i/3+1)%L_ == 0)
            temp[2] -= 3*L_;
        if ( ( (i/3)%(L_*L_*L_) >= (L_*L_*L_) - L_*L_) && ( (i/3)%(L_*L_*L_) < (L_*L_*L_) ) )
            temp[3] -= 3*L_*L_*L_;

        Plaquette.push_back(temp);

        //YZ plane
        temp[0] = i+1;
        temp[1] = i+2;
        temp[2] = i+3*L_ + 2;
        temp[3] = i+3*L_*L_ + 1;
        if ( ( (i/3)%(L_*L_) >= (L_*L_) - L_) && ( (i/3)%(L_*L_) < (L_*L_) ) )
            temp[2] -= 3*L_*L_;
        if ( ( (i/3)%(L_*L_*L_) >= (L_*L_*L_) - L_*L_) && ( (i/3)%(L_*L_*L_) < (L_*L_*L_) ) )
            temp[3] -= 3*L_*L_*L_;

        Plaquette.push_back(temp);

    }//i

    //DEBUG: check if Plaquette has any errors
    vector<int> Check(Plaquette.size(),0);
    //cout<<"Check size : "<<Check.size()<<endl;
    for (int j=0; j<Check.size(); j++)
        for (int k=0; k<Plaquette[j].size(); k++)
            Check[Plaquette[j][k]]++;

    for (int j=0; j<Check.size(); j++)
        if (Check[j] != 4) cout<<"Plaquette error \n";
    //cout<<j<<" "<<Check[j]<<endl;

    Energy = CalcEnergy(sigma);      
    cout<<Energy<<endl;      

    //Now, make the data structure used to relate the DOF to the 4 plaquettes
    All_Neighbors.resize(N1);
    for (int i=0; i<Plaquette.size(); i++)
        for (int j=0; j<Plaquette[i].size(); j++)
            All_Neighbors[Plaquette[i][j]].push_back(i);


}//constructor


//print
void ThreeD12Code::print(){

    cout<<L_<<" "<<D_<<" "<<N1<<endl;

    cout<<"Plaquette \n";
    for (int i=0; i<Plaquette.size(); i++){
        //cout<<i<<" ";
        for (int j=0; j<4; j++)
            cout<<Plaquette[i][j]<<" ";
        cout<<endl;
    }//i

    for (int i=0; i<All_Neighbors.size(); i++){
        for (int j=0; j<All_Neighbors[i].size(); j++)
            cout<<All_Neighbors[i][j]<<" ";
        cout<<endl;
    }

}//print


//loops through to calculate the energy
double ThreeD12Code::CalcEnergy(Spins & sigma){

    double eTemp = 0.0;

    for (int i=0; i<Plaquette.size(); i++){
        eTemp -=  sigma.spin[Plaquette[i][0]]*sigma.spin[Plaquette[i][1]]
            *sigma.spin[Plaquette[i][2]]*sigma.spin[Plaquette[i][3]];
    }//i

    return eTemp;

}

//the fast way to calculte the new energy
double ThreeD12Code::CalcEnergyDiff(Spins & sigma, const int & flipsite){

    double DeltaE = 0.0;
    double spinProd;

    for (int j=0; j<All_Neighbors[flipsite].size(); j++){
        spinProd = 1; 
        for(int k=0; k<Plaquette[0].size(); k++) {
            spinProd *= sigma.spin[ Plaquette[All_Neighbors[flipsite][j]][k] ];
        }//k

        DeltaE += -spinProd; //ferromagnetic
    }//j

    DeltaE *= 2.0; //double counting

    return DeltaE;

}

//Calculates a number of single-spin flips
void ThreeD12Code::LocalUpdate(Spins & sigma, const double & T, MTRand & ran){

    int site;  //random site for update
    double Eold, Enew, Ediff;
    double m_rand; //metropolis random number

    for (int j=0; j<N1; j++){ //peform N random single spin flips

        site = ran.randInt(N1-1);
        //cout<<"site is "<<site<<endl;

        sigma.flip(site);  //trial flip
        Eold = Energy;
        //Enew = CalcEnergy(sigma); //slow way
        //Ediff = Enew - Eold;
        Ediff = CalcEnergyDiff(sigma,site); //fast way
        Enew = Eold + Ediff;

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


