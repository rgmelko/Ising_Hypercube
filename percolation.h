#ifndef PERCOLATION_H
#define PERCOLATION_H

// percolation.h: a class that performs percolation measurements given
//   1) The arbitrary list of Neighbors or Nodes
//   2) The "occupation" of the nodes, as 0 or 1
//
// Hoshen-Kopelman algorithm by Grant Watson, see https://github.com/MelkoCollective/extended_hk
// last update August 17, 2013, 5892ff451f1f014b6c6bc76681476a2e84527d51

#include <boost/multi_array.hpp>

#define PRINT_RED(x) std::cout << "\033[1;31m" << x << "\033[0m" << " "
#define PRINT_BLUE(x) std::cout << "\033[1;34m" << x << "\033[0m" << " "
#define PRINT_GREEN(x) std::cout << "\033[1;32m" << x << "\033[0m" << " "
#define PRINT_YELLOW(x) std::cout << "\033[1;33m" << x << "\033[0m" << " "

class Percolation
{

    private:
        int N_;
        boost::multi_array<int, 1> UniqueClusters;
        boost::multi_array<int, 1> ClustSize;  //the size of each cluster

        //all below from Grant's hk.cpp -------------
        int n_labels; //the length of the labels array
        int *labels;
        int uf_find(int x);
        int uf_union(int x, int y);
        int uf_make_set(void);
        void uf_initialize(int max_labels);
        void uf_done(void);
        //-------------------------------------------

    public:
        double Avg_Clust_Size;   
        double Largest_Clust_Size;   
        int NumberClusters;

        Percolation(const int & N);    
        void zero();
        void print();
        void DetermineClusters(const boost::multi_array<int, 2>& nbs, 
                const boost::multi_array<int, 1>& occupancy);
        //void record(double & energy, Spins & sigma);
        void output(const double & T, const int & MCS);

        //Headers from Grant's hk.h file, worked into this class
        //-------------------------------------------
        void extended_hoshen_kopelman(boost::multi_array<int, 1>& node_labels,
                const boost::multi_array<int, 2>& nbs,
                const boost::multi_array<int, 1>& occupancy);

        void extended_hk_no_boost(int* node_labels, int const* const* nbs,
                const int* occupancy, int N, int m);

        int hoshen_kopelman(int **matrix, int m, int n);

};

//constructor
Percolation::Percolation(const int & N){

    N_ = N; //number of ising variables that can be in a cluster

    n_labels = 0; /* length of the labels array */

    Avg_Clust_Size= 0;
    Largest_Clust_Size = 0;

}

void Percolation::print(){

    PRINT_YELLOW("Cluster determination")<<endl;
    cout<<"there are "<<NumberClusters<<" labels \n";
    for (int i=0; i<N_; i++){
        cout<<UniqueClusters[i]<<" ";
    }
    cout<<endl;

    //print each cluster's size
    PRINT_GREEN("Each cluster has a size (zero included):")<<endl;
    for (int i=0; i<=NumberClusters; i++)
        cout<<ClustSize[i]<<endl;

    cout<<Avg_Clust_Size<<" ";
    cout<<Largest_Clust_Size<<endl;

}//print

void Percolation::zero(){

    Avg_Clust_Size = 0.0;
    Largest_Clust_Size= 0.0;

}//zero

void Percolation::DetermineClusters(const boost::multi_array<int, 2>& nbs, 
        const boost::multi_array<int, 1>& occupancy) {

    extended_hoshen_kopelman(UniqueClusters,nbs,occupancy);

    ClustSize.resize(boost::extents[NumberClusters+1]); //resize
    fill(ClustSize.begin(),ClustSize.end(),0); //reinitialize

    for (int i=0; i<N_; i++) { //note 0 is not a cluster
        ClustSize[UniqueClusters[i]] ++;
    }

    double Asize=0;
    int Lclust =0;
    if (NumberClusters !=0 ){
        for (int i=1; i<(NumberClusters+1); i++){
            Asize += 1.0*ClustSize[i]; //average cluster size
            if (ClustSize[i] > Lclust) Lclust = ClustSize[i]; //largest cluster size
        }
        Asize /= 1.0*NumberClusters;
    }

    Avg_Clust_Size += Asize;
    Largest_Clust_Size += 1.0*Lclust;


}//DetermineClusters

void Percolation::output(const double & T, const int & MCS){

    ofstream cfout;
    cfout.open("01.data",ios::app);

    cfout<<T<<" ";
    cfout<<Avg_Clust_Size/(1.0*MCS*N_)<<" ";
    cfout<<Largest_Clust_Size/(1.0*MCS*N_)<<" ";
    cfout<<endl;

    cfout.close();

}//output




/*-------------------------------------------------------------------------------------------
  below, the hk.cpp file from 
https://github.com/MelkoCollective/extended_hk
-------------------------------------------------------------------------------------------*/

/* Tobin Fricke's original code with an extended version of the HK algorithm
   added. The extension borrows some ideas from:
http://gaia.pge.utexas.edu/papers/AFTWPPhysicaA.pdf

Requirements:
-Boost
-compiler ? - ?

Copyright (c) September 9, 2000, by Tobin Fricke <tobin@pas.rochester.edu>

Extended 2013-08-14 Grant Watson
Modified 2002-03-09 Tobin Fricke
Modified substantially 2004-04-21 by Tobin Fricke

Fricke's code available at:
http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
 */

#include <boost/phoenix.hpp>
#include <algorithm>
#include <cassert>
#include <new>

using namespace std;

/* Implementation of Union-Find Algorithm */

/* The 'labels' array has the meaning that labels[x] is an alias for the label x; by
   following this chain until x == labels[x], you can find the canonical name of an
   equivalence class.  The labels start at one; labels[0] is a special value indicating
   the highest label already used. */

// SEE PRIVATE CLASS MEMBERS
//int *labels;
//int n_labels = 0; /* length of the labels array */

/*  uf_find returns the canonical label for the equivalence class containing x */

int Percolation::uf_find(int x) {
    int y = x;
    while (labels[y] != y)
        y = labels[y];

    while (labels[x] != x) {
        int z = labels[x];
        labels[x] = y;
        x = z;
    }
    return y;
}

/*  uf_union joins two equivalence classes and returns the canonical label of the resulting class. */

int Percolation::uf_union(int x, int y) {
    return labels[uf_find(x)] = uf_find(y);
}

/*  uf_make_set creates a new equivalence class and returns its label */

int Percolation::uf_make_set(void) {
    labels[0]++;
    assert(labels[0] < n_labels);
    labels[labels[0]] = labels[0];
    return labels[0];
}

/*  uf_intitialize sets up the data structures needed by the union-find implementation. */

void Percolation::uf_initialize(int max_labels) {
    n_labels = max_labels;
    labels = new int[n_labels];
    labels[0] = 0;
}

/*  uf_done frees the memory used by the union-find data structures */

void Percolation::uf_done(void) {
    n_labels = 0;
    delete[] labels;
    labels = 0;
}

/* End Union-Find implementation */

/* A generalized version of the HK algorithm for arbitrary networks of nodes.
 *
 * INPUT:
 * -nbs: 2D array. ith row is the neighbours of node i. If a node has
 *       less neighbours than # of cols, fill the rest of the row with -1.
 * -occupancy: 1D array with the occupation number (0 or 1) of the nodes.
 * OUTPUT:
 * -node_labels: an array of the labels of the nodes.
 */
void Percolation::extended_hoshen_kopelman(boost::multi_array<int, 1>& node_labels,
        const boost::multi_array<int, 2>& nbs,
        const boost::multi_array<int, 1>& occupancy) {
    // Typedefs
    typedef boost::multi_array<int, 2> array_2t;
    typedef boost::multi_array<int, 1> array_1t;
    typedef boost::multi_array_types::index_range range_t;
    typedef boost::multi_array<int, 1>::const_iterator c_iter;

    // Tools for multi-arrays.
    array_1t::extent_gen extents;
    array_2t::index_gen indices;

    // Number of nodes.
    const int N = nbs.shape()[0];
    const int m = nbs.shape()[1];

    // Initialize node_labels with N+1 since labels live in [1,N].
    const int unlabelled = N+1;
    node_labels.resize(extents[N]);
    fill(node_labels.begin(),node_labels.end(),N+1);

    // Initialize memory for binary forest of labels.
    uf_initialize(N);

    // Iterate over nodes and perform clustering.
    for (int i = 0; i < N; ++i) {
        if (occupancy[i]) {

            // Get neighbours of node i ('i'th row of neighbours)
            int n_nbs = m;
            while (nbs[i][n_nbs - 1] == -1) // Find last neighbour.
                n_nbs--;
            array_2t::const_array_view<1>::type node_nbs =
                nbs[ indices[i][range_t(0,n_nbs)] ];

            // Get subset of labels using node_nbs as indices.
            array_1t node_nbs_labels(extents[n_nbs]);
            for (unsigned int j = 0; j < n_nbs; ++j) {
                node_nbs_labels[j] = node_labels[node_nbs[j]];
            }
            c_iter begin = node_nbs_labels.begin();
            c_iter end = node_nbs_labels.end();

            // Check if node has no labelled neighbours.
            bool is_alone = 1;
            for (c_iter it = begin; is_alone && (it != end); ++it){
                is_alone *= (*it == unlabelled);
            }

            // Labelling + merging
            if(is_alone)
                node_labels[i] = uf_make_set();
            else {
                // Find smallest label of the neighbours.
                int min_label = *min_element(begin, end);

                // Apply the minimum label to all labelled neighbours + current node.
                node_labels[i] = min_label;
                for (c_iter nb = node_nbs.begin(); nb != node_nbs.end(); ++nb)
                    if (node_labels[*nb] != unlabelled)
                        uf_union(min_label, node_labels[*nb]);
            }

        } //occupancy
    } //node

    /* This is a little bit sneaky.. we create a mapping from the canonical labels
       determined by union/find into a new set of canonical labels, which are
       guaranteed to be sequential. */

    NumberClusters = 0;
    int *new_labels = new int[n_labels]; // allocate array, initialized to zero
    for (int i = 0; i < n_labels; i++) new_labels[i] = 0;

    for (int i = 0; i < N; i++)
        if (occupancy[i]) {
            int x = uf_find(node_labels[i]);
            if (new_labels[x] == 0) {
                new_labels[0]++;
                new_labels[x] = new_labels[0];
                NumberClusters ++;
            }
            node_labels[i] = new_labels[x];
        }
        else {
            node_labels[i] = 0; // Replace placeholders with 0.
        }

    // Cleanup
    delete[] new_labels;
    uf_done();
}

#endif
