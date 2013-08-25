#ifndef MEASURE_H
#define MEASURE_H

// measure.h: a class that performs statistical measurements of estimators
// Energy, Specific Heat, Magnetization, Susceptibility

#include "spins.h"

typedef boost::multi_array<int, 2> array_2t;

class Measure
{
    private:
       int Nspin;
       int MCS;

    public:
      double TOT_energy;   //energy
      double TOT_energy2;  //energy^2
      double TOT_Mag;    //magnetization s
      double TOT_Mag2;    //magnetization squared
      double TOT_WilX;    //magnetization squared

      //vector <double> Wilson_RunningAvg;

      Measure(const int &, const PARAMS &);
      void zero();
      void record(double & energy, Spins & sigma, const array_2t &);
      void output(const double &);
      void outputWilsonLoop(const Spins & , const array_2t & , const int & );
 
};

//constructor
Measure::Measure(const int & N, const PARAMS & p){

    Nspin = N;
    MCS = p.MCS_;

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_Mag = 0.0;
    TOT_Mag2 = 0.0;
    TOT_WilX= 0.0;

	//Wilson_RunningAvg.assign(p.Dim_,0.0);
}

//zero
void Measure::zero(){

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_Mag = 0.0;
    TOT_Mag2 = 0.0;
    TOT_WilX= 0.0;	

	//for (int d=0; d<Wilson_RunningAvg.size(); d++)
    //Wilson_RunningAvg[d] = 0.0;

}


void Measure::record(double & energy, Spins & sigma, const array_2t & WilsonLoops){

    TOT_energy += energy;
    TOT_energy2 += energy * energy;

    long int mag = 0;
    for (int i=0; i<sigma.spin.size(); i++)
        mag += sigma.spin[i];

    TOT_Mag += 1.0*mag;
    TOT_Mag2 += 1.0*mag*mag;

    int prod = 1;
	int L = WilsonLoops[0].size();
    for (int i=0; i<L; i++)
		prod *= sigma.spin[WilsonLoops[0][i]];

    TOT_WilX += 1.0*prod;

}//update

void Measure::output(const double & T){

	ofstream cfout;
	cfout.open("00.data",ios::app);

    cfout<<T<<" ";
    cfout<<TOT_energy/(1.0*MCS * Nspin)<<" ";
	double Cv = TOT_energy2/(1.0*MCS) - TOT_energy*TOT_energy/(1.0*MCS*MCS); 
    cfout<<Cv/(T*T*1.0*Nspin)<<" ";
    cfout<<TOT_Mag2/(1.0*MCS * Nspin*Nspin)<<" ";
	double susc = TOT_Mag2/(1.0*MCS) - TOT_Mag*TOT_Mag/(1.0*MCS*MCS); 
    cfout<<susc/(T*1.0*Nspin)<<" ";
    cfout<<TOT_WilX/(1.0*MCS)<<"\n";

	cfout.close();

}//output

void Measure::outputWilsonLoop(const Spins & sigma, const array_2t & WilsonLoops, const int & SimNum){

    int SimTemp = SimNum;

    char fname[8];

    if (SimTemp/10 == 0)
        fname[0]='0';
    else if (SimTemp/20 == 0) {
        fname[0]='1';
        SimTemp = SimTemp - 10;   }
    else if (SimTemp/30 == 0) {
        fname[0]='2';
        SimTemp = SimTemp - 20;   }
    else if (SimTemp/40 == 0) {
        fname[0]='3';
        SimTemp = SimTemp - 30;   }

    if (SimTemp == 0) fname[1] = '0';
    else if (SimTemp%9 == 0) fname[1] = '9';
    else if (SimTemp%8 == 0) fname[1] = '8';
    else if (SimTemp%7 == 0) fname[1] = '7';
    else if (SimTemp%6 == 0) fname[1] = '6';
    else if (SimTemp%5 == 0) fname[1] = '5';
    else if (SimTemp%4 == 0) fname[1] = '4';
    else if (SimTemp%3 == 0) fname[1] = '3';
    else if (SimTemp%2 == 0) fname[1] = '2';
    else if (SimTemp%1 == 0) fname[1] = '1';

    fname[2] = '.';
    fname[3] = 'w';
    fname[4] = 'n';
    fname[5] = 'u';
    fname[6] = 'm';
    fname[7] = '\0';


	int L = WilsonLoops[0].size();
	int D = WilsonLoops.size();

 	ofstream cfout;
	cfout.open(fname,ios::app);
   
    double prod;
	for (int d=0; d<D; d++){
		prod = 1;
		for (int i=0; i<L; i++)
			prod *= sigma.spin[WilsonLoops[d][i]];
		//cfout<<prod<<" ";
		//Wilson_RunningAvg[d] += prod;
       //cfout<<Wilson_RunningAvg[d]/(1.0*MCstep)<<" ";
       cfout<<prod<<" ";
	}
	cfout<<endl;

	cfout.close();

}//outputWilsonLoop


#endif
