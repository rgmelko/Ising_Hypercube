#ifndef MEASURE_H
#define MEASURE_H

#include "spins.h"

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

      Measure(const int &, const PARAMS &);
      void zero();
      void record(double & energy, Spins & sigma);
      void output(const double &);
  
};

//constructor
Measure::Measure(const int & N, const PARAMS & p){

    Nspin = N;
    MCS = p.MCS_;

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_Mag = 0.0;
    TOT_Mag2 = 0.0;
}

//zero
void Measure::zero(){

    TOT_energy = 0.0;
    TOT_energy2 = 0.0;
    TOT_Mag = 0.0;
    TOT_Mag2 = 0.0;
}


void Measure::record(double & energy, Spins & sigma){

    TOT_energy += energy;
    TOT_energy2 += energy * energy;

    long int mag = 0;
    for (int i=0; i<sigma.spin.size(); i++)
        mag += sigma.spin[i];

    TOT_Mag += 1.0*mag;
    TOT_Mag2 += 1.0*mag*mag;


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
    cfout<<susc/(T*1.0*Nspin)<<"\n";

	cfout.close();

}//output


#endif
