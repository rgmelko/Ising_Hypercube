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
        int N_;   //number of lattice sites
        int D_;   //Dimension
        int L_;   //Linear size

        double Energy;  //total energy of the system

        //All the 1-cells (bonds) that are attached to 2-cells (faces)
        //vector<vector<int> > All_Neighbors; 

        //The Face operators
        vector<vector<int> > Plaquette;
        ThreeD12Code(Spins & sigma, HyperCube & cube); 
        double CalcEnergy(Spins & sigma);
        void LocalUpdate(Spins & sigma, double & T, MTRand & ran);

        void print();

    private:

        int Faces; //The total number of faces

};

//constructor
ThreeD12Code::ThreeD12Code(Spins & sigma, HyperCube & cube){

    L_ = cube.L_;
    D_ = 3;
    N_ = 3*cube.N_; //3D lattice BOND variables...

	sigma.resize(N_); //these are the degrees of freedom

    Faces = N_;  //in 3D the number of faces equals the numbers of DOFs

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
    }//i

    //cout<<CalcEnergy(sigma)<<endl;      

}//constructor


//print
void ThreeD12Code::print(){

    cout<<L_<<" "<<D_<<" "<<N_<<endl;

//    for (int i=0; i<All_Neighbors.size(); i++){
//        cout<<i<<" ";
//        for (int j=0; j<Bonds_Per_Site; j++)
//            cout<<All_Neighbors[i][j]<<" ";
//        cout<<endl;
//    }//i

    cout<<"Plaquette \n";
    for (int i=0; i<Plaquette.size(); i++){
        cout<<i<<" ";
        for (int j=0; j<4; j++)
            cout<<Plaquette[i][j]<<" ";
        cout<<endl;
    }//i

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

//Calculates a number of single-spin flips
void ThreeD12Code::LocalUpdate(Spins & sigma, double & T, MTRand & ran){

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


