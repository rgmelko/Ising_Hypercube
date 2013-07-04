#ifndef HYPERCUBE_H
#define HYPERCUBE_H

// hypercube.h
// a class to make a hypercubic lattice of arbitrary dimension

#define PRINT_RED(x) std::cout << "\033[1;31m" << x << "\033[0m" << " "
#define PRINT_BLUE(x) std::cout << "\033[1;34m" << x << "\033[0m" << " "
#define PRINT_GREEN(x) std::cout << "\033[1;32m" << x << "\033[0m" << " "
#define PRINT_YELLOW(x) std::cout << "\033[1;33m" << x << "\033[0m" << " "

#include <vector>
#include <iostream>

using namespace std;

class HyperCube
{
    public:
        int L_; //linear size
        int D_; //dimension
        int N_; //total number of sites

        //the lattice is a vector of vectors: no double counting
        vector<vector<int> > Neighbors;
        vector<vector<int> > Coordinates;

        //public functions
        HyperCube(int L, int D);
        void print();
         
    private:
        int myPow(int, int);

};

//constructor
//takes the lattice linear size L, and dimension D
HyperCube::HyperCube(int L, int D){

    L_ = L;
    D_ = D;
    N_ = myPow(L_,D_);

    //build the nearest-neighbor connections
    vector<int> temp;  //the "inner" vector
    int pair; //the bond pair index
    for (int i=0;i<N_;i++){

        temp.clear();
        for (int j=0;j<D_;j++){
            if ( ( (i+myPow(L,j))%(myPow(L,(j+1))) >= 0) 
                    && ( (i+myPow(L,j))%(myPow(L,(j+1))) < myPow(L,j) ) )
                pair = i+myPow(L,j) - myPow(L,(j+1));
            else
                pair = i+myPow(L,j);

            temp.push_back(pair);
        }

        Neighbors.push_back(temp);
    }//i

    //build the (x,y,z,...) coordinates of each lattice site
    temp.clear();
    temp.assign(D_,0);  //D integers with value 0

    for (int i=0;i<N_;i++){
        Coordinates.push_back(temp);

        if ( (temp[0]+1) % L_ == 0){ //end of x-row

            temp[0] = 0; //reset

            for (int j=1;j<D_;j++){
                if ( (temp[j]+1) % L_ == 0)
                    temp[j] = 0; //reset
                else{
                    temp[j]++;
                    break;
                }
            }//j
        }//if
        else
            temp[0]++;
    }


}//constructor

//a print function
void HyperCube::print(){

    cout<<"L D N \n";
    cout<<L_<<" "<<D_<<" "<<N_<<endl;

    cout<<"Neighbor list:"<<endl;
    for (int i=0;i<Neighbors.size();i++){
        cout<<i<<" ";
        for (int j=0;j<Neighbors[i].size();j++){
            //cout<<j<<" ";
            //cout<<Neighbors[i][j]<<" ";
            PRINT_RED(Neighbors[i][j]);
        }
        cout<<endl;
    }//i

    cout<<"Coordinates:"<<endl;
    for (int i=0;i<Coordinates.size();i++){
        cout<<i<<" ";
        for (int j=0;j<Coordinates[i].size();j++){
            //cout<<j<<" ";
            //cout<<Coordinates[i][j]<<" ";
			PRINT_YELLOW(Coordinates[i][j]);

        }
        cout<<endl;
    }//i

}

//a simple function to calculate powers of an integer: not for general use
int HyperCube::myPow (int x, int p) {
  int i = 1;
  for (int j = 1; j <= p; j++)  i *= x;
  return i;
}

#endif