#ifndef GENERALD_1_2_CODE_H
#define GENERALD_1_2_CODE_H

// fourD_1_2_code.h
// a class to perform a simple metropolis MC on a 2D (1,2) toric code

#define PRINT_RED(x) std::cout << "\033[1;31m" << x << "\033[0m" << " "
#define PRINT_BLUE(x) std::cout << "\033[1;34m" << x << "\033[0m" << " "
#define PRINT_GREEN(x) std::cout << "\033[1;32m" << x << "\033[0m" << " "
#define PRINT_YELLOW(x) std::cout << "\033[1;33m" << x << "\033[0m" << " "

#include "spins.h"
#include "MersenneTwister.h"
#include <vector>
#include <iostream>

using namespace std;

typedef boost::multi_array<int, 2> array_2t;
typedef boost::multi_array<int, 1> array_1t;

class GeneralD12Code
{
    private:
        array_2t dims2plane;
        int Nplane;   //number of planes (Nchoose2)        

    public:
        int N0;
        int N1;   //number of DEGREES OF FREEDOM
        int N2;   //number of 2 cells
        int N3;
        int D_;   //Dimension
        int L_;   //Linear size

        double Energy;  //total energy of the system

        //All of the 2-cells which are attached to a given 1-cell
        vector<vector<int> > All_Neighbors; 

        //The neighbor list for 2-cells: defined if sharing a 3-cell
        array_2t TwoCellNeighbors;

        //the occupancy - for percolation
        array_1t occupancy;

        //The Face operators
        vector<vector<int> > Plaquette;
        vector<vector<int> > Cubes;
        GeneralD12Code(Spins & sigma, HyperCube & cube); 
        double CalcEnergy(Spins & sigma);
        double CalcEnergyDiff(Spins & sigma, const int & flipsite);
        void CalculateOccupancy(const Spins & sigma);
        void PreparePercolation(const Spins & sigma, const HyperCube & cube);

        void LocalUpdate(Spins & sigma, const double & T, MTRand & ran);

        void print();


};

//constructor
GeneralD12Code::GeneralD12Code(Spins & sigma, HyperCube & cube){

    L_ = cube.L_;
    D_ = cube.D_;
    N0 = cube.N_;
    N1 = D_*cube.N_; //number of 1 cells

    sigma.resize(N1); //these are the degrees of freedom (1 cells)
    sigma.randomize();

    //use it to built the sigma-z plaquettes
    vector <int> temp;
    temp.assign(4,0);  //assign 4 zeros to this vector

    //ANN's IDEA JUST SO YOU KNOW
    //Input: plaquette#
    //Output: the 4 1-cells associated with that plaquette
    for (int v=0; v<N0; v++ ){ //loop over 0-cells

        for (int i=0; i<(D_-1); i++){ //loop that defines all 2-cells per vertex
            for (int j=0; j<D_; j++){

                if (i<j){  // cout<<i<<" "<<j<<endl;

                    temp[0] = D_*v+i;
                    temp[1] = D_*v+j;
                    temp[2] = D_*cube.Neighbors[v][i]+j;
                    temp[3] = D_*cube.Neighbors[v][j]+i;
                    Plaquette.push_back(temp);

                }//if

            }//j
        }//i

    }//v

    N2 = Plaquette.size(); //number of 2 cells

    //DEBUG: check if Plaquette has any errors
    //vector<int> Check(Plaquette.size(),0);
    vector<int> Check(N1,0);
    //cout<<"Check size : "<<Check.size()<<endl;
    for (int j=0; j<Plaquette.size(); j++)
        for (int k=0; k<Plaquette[j].size(); k++)
            Check[Plaquette[j][k]]++;

    for (int j=0; j<Check.size(); j++){
        if (Check[j] != 2*(D_-1)){ 
            cout<<"Plaquette error \n";
            cout<<j<<" "<<Check[j]<<endl;
        }   
    }

    //Now, make the data structure used to relate the DOF to the 4 plaquettes
    All_Neighbors.resize(N1);
    for (int i=0; i<Plaquette.size(); i++)
        for (int j=0; j<Plaquette[i].size(); j++)
            All_Neighbors[Plaquette[i][j]].push_back(i);


    Energy = CalcEnergy(sigma);      
    cout<<"Energy: "<<Energy<<endl;      


}//constructor


//Prepare data structures for percolation: valid in D>2 only
void GeneralD12Code::PreparePercolation(const Spins & sigma, const HyperCube & cube){

    if (D_ < 3) cout<<"ERROR: DIMENSION TOO LOW FOR PERCOLATION \n";

    //ANN's OTHER IDEA (JUST SO YOU KNOW)
    //Create an object that translates between 2 dimensions
    //and the plane they represent
    dims2plane.resize(boost::extents[D_][D_]);
    int tempCount=0;
    for(int i=0; i<D_; i++){
        for(int j=0; j<D_; j++){
            if(i<j){ 
                dims2plane[i][j]=tempCount;
                dims2plane[j][i]=tempCount;
                tempCount++;
            }
            else{ dims2plane[i][j]=-99; }
        }//j
    }//i
    Nplane = tempCount; //this is the number of diff planes (Nchoose2)
    //----Done filling dims2plane

    //Now creating Cubes
    //Input: cube identifier #
    //Output: the 6 2-cells associated with that cube
    vector <int> temp;
    temp.assign(6,0);
    for (int v=0; v<N0; v++ ){ //loop over 0-cells

        for (int i=0; i<D_; i++){ //loop that defines all 2-cells per vertex
            for (int j=0; j<D_; j++){
                for (int k=0; k<D_; k++){

                    if ((i<j)&&(j<k)){   //cout<<i<<" "<<j<<" "<<k<<endl;

                        temp[0] = Nplane*v + dims2plane[i][j];
                        temp[1] = Nplane*v + dims2plane[i][k];
                        temp[2] = Nplane*v + dims2plane[j][k];
                        temp[3] = Nplane*cube.Neighbors[v][k] + dims2plane[i][j];
                        temp[4] = Nplane*cube.Neighbors[v][j] + dims2plane[i][k];
                        temp[5] = Nplane*cube.Neighbors[v][i] + dims2plane[j][k];
                        Cubes.push_back(temp);

                    }//if
                }//k
            }//j
        }//i
    }//v

    N3 = Cubes.size();//Number of 3-cells
    cout << "N3=" << N3 <<endl;

    occupancy.resize(boost::extents[N2]); //calculate percolation objects

    CalculateOccupancy(sigma); //for percolation

    //Below defines which 2-cells are neighbors: belong to the same 3-cell (for percolation)


    TwoCellNeighbors.resize(boost::extents[N2][10*(D_-2)]); 
    //initialize
    for (int i=0; i<TwoCellNeighbors.size(); i++)
        for (int j=0; j<TwoCellNeighbors[i].size(); j++)
            TwoCellNeighbors[i][j] = -99;

    int p1, p2;
    for (int v3=0; v3<Cubes.size(); v3++){

        for (int i=0; i<Cubes[v3].size(); i++){
            p1 = Cubes[v3][i];
            for (int j=0; j<Cubes[v3].size(); j++){
                p2 = Cubes[v3][j];
                if (p1 != p2){

                    for (int k=0; k<TwoCellNeighbors[p1].size(); k++) //no push_back
                        if (TwoCellNeighbors[p1][k] == -99){
                            TwoCellNeighbors[p1][k] = p2;
                            break;
                        }

                }//if
            }//j
        }//i

    }//v3

}//PreparePercolation


//print
void GeneralD12Code::print(){

    cout<<L_<<" "<<D_<<" "<<N1<<" "<<N2<<endl;

    cout<<"Plaquette \n";
    for (int i=0; i<Plaquette.size(); i++){
        PRINT_RED(i);
        for (int j=0; j<4; j++)
            cout<<Plaquette[i][j]<<" ";
        //PRINT_RED(Plaquette[i][j]);
        cout<<endl;
    }//i

    for (int i=0; i<All_Neighbors.size(); i++){
        PRINT_GREEN(i);
        for (int j=0; j<All_Neighbors[i].size(); j++){
            cout<<All_Neighbors[i][j]<<" ";
        }
        //PRINT_GREEN(All_Neighbors[i][j]);
        cout<<endl;
    }

//    for (int i=0; i<N3; i++){
//        //PRINT_BLUE(i);
//        for (int j=0; j<6; j++){
//            cout<<Cubes[i][j]<<" ";
//        }
//        cout<<endl;
//    }//i
//
//    cout<<endl;
//
//    for (int i=0; i<N2; i++){
//        //PRINT_RED(i);
//        for (int j=0; j<10*(D_-2); j++){
//            cout<<TwoCellNeighbors[i][j]<<" ";
//        }
//        cout<<endl;
//    }//i


}//print


//loops through to calculate the energy
double GeneralD12Code::CalcEnergy(Spins & sigma){

    double eTemp = 0.0;

    for (int i=0; i<Plaquette.size(); i++){
        eTemp -=  sigma.spin[Plaquette[i][0]]*sigma.spin[Plaquette[i][1]]
            *sigma.spin[Plaquette[i][2]]*sigma.spin[Plaquette[i][3]];
    }//i

    return eTemp;

}

//loops through to calculate the occupancy for percolation
void GeneralD12Code::CalculateOccupancy(const Spins & sigma){

    double eTemp = 0.0;

    int no_defect;
    for (int i=0; i<Plaquette.size(); i++){

        no_defect = sigma.spin[Plaquette[i][0]]*sigma.spin[Plaquette[i][1]]
            *sigma.spin[Plaquette[i][2]]*sigma.spin[Plaquette[i][3]];

        if (no_defect == -1) occupancy[i] = 1;  //this is a defect
        else occupancy[i] = 0;

        eTemp -= 1.0*no_defect;

    }//i

    if (eTemp != Energy) cout<<"Plaquette Energy Problem  \n";

}//CalculateOccupancy


//the fast way to calculte the new energy
double GeneralD12Code::CalcEnergyDiff(Spins & sigma, const int & flipsite){

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
void GeneralD12Code::LocalUpdate(Spins & sigma, const double & T, MTRand & ran){

    int site;  //random site for update
    double Eold, Enew, Ediff;
    double m_rand; //metropolis random number

    bool accept;
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
        accept = true;
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
                accept = false;
            }
        }//Metropolis

        if (accept == true) { //update the occupancy list since the flip was accepted
            for (int k=0; k<All_Neighbors[site].size(); k++) //bit flip
                occupancy[All_Neighbors[site][k]] = 1 -occupancy[All_Neighbors[site][k]]; 
        }//accept

    }//j

    //cout<<"Emod "<<Energy<<endl;
}//LocalUpdate



#endif


