/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){

  int M = 250000;     // realizations
  int N = 50;         // blocks
  string filename = "../../Results/gofr_NVE_liquid.dat";

  Input();                                          // input
  if(restart=="true") Equilibrate_system(3000);     // equilibration if needed
  blocking_on_MD(M, N, filename);                   // blocking
  ConfFinal();

  return 0;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf, ReadOld;
  ofstream SaveConf;
  double ep, ek, pr, et, vir;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;
  ReadInput >> npart;
  ReadInput >> rho;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;

  if(restart=="false"){
    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
    cout << "Number of particles = " << npart << endl;
    cout << "Density of particles = " << rho << endl;
    cout << "Volume of the simulation box = " << vol << endl;
    cout << "Edge of the simulation box = " << box << endl;
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
  }
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

// for evaluating the g(r)

  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;     // updating cardinality of my measure
  bin_size = (box/2.0)/(double)nbins;

  stima_gofr.resize(nbins);

//Read r(t) positions

  if(restart=="true"){
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

    SaveConf.open("old.0", ios::out | ios::trunc);
    //Prepare initial velocities
    cout << "Restart is True: preparing first time configuration in the new old.0 file" << endl;
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = rand()/double(RAND_MAX) - 0.5;
      vy[i] = rand()/double(RAND_MAX) - 0.5;
      vz[i] = rand()/double(RAND_MAX) - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
      SaveConf << xold[i]/box << "  " << yold[i]/box << "  " << zold[i]/box << endl;
    }
    SaveConf.close();

  } else if (restart=="false"){
      cout << "Restart is false: opening old.0 and using it for defining r(t-dt)" << endl;
      ReadOld.open("old.final");
      ReadConf.open("config.final");

      for (int i=0; i<npart; ++i){
        ReadOld  >> xold[i] >> yold[i] >> zold[i];
        ReadConf >> x[i] >> y[i] >> z[i];
        xold[i] = xold[i] * box;
        yold[i] = yold[i] * box;
        zold[i] = zold[i] * box;
        x[i] = x[i]*box;
        y[i] = y[i]*box;
        z[i] = z[i]*box;
      }
      ReadConf.close();
      ReadOld.close();
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ THERMALIZATION PHASE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Equilibrate_system(int N){

  // SAME COMMENTS DONE IN LAB_04

  cout << "####################################################################" << endl;
  cout << "Thermalization phase of the simulation." << endl;
  cout << "Running "<< N << " steps that will be ignored at the end of this phase." << endl << endl;
  cout << "####################################################################" << endl;

  int hist_dimension = 100;
  int hist_count     = 0;
  vector<double> mean_v2_history;

  for(int i=0; i<N; i++){

    Move();
    if(hist_count > 220) mean_v2_history.push_back(eval_mean_v2());

    //if(i==int(N/10)){temp = 1.0; cout << "Changing temperature T*: 0.8 --> 1.0" << endl;}
    //if(i==int(N/5)) {temp = 1.2; cout << "Changing temperature T*: 1.0 --> 1.2" << endl;}

    if((i+1)%(300)==0){
      cout << "Thermalization process is running, step " << i+1 << "/" << N << ". Rescaling velocities." << endl;
      cout << "The kinetic energy is evaluated using the last" << mean_v2_history.size() << " istant measures." << endl;
      rescale_velocities(mean_v2_history);
      mean_v2_history.clear();
      hist_count = 0;
    }
    if(i%100==0) Measure(true);
    hist_count++;

  }
  set_restart("true","false");
  cout << "This is the end of the thermalization phase, let's the simulation begin. " << endl;
  cout << "####################################################################" << endl << endl;

  ConfFinal();

  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FZ FOR REPLACING RESTART VALUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~

void set_restart(string current, string new_val){
  string strReplace = current;
  string strNew = new_val;
  ifstream filein("input.dat");
  ofstream fileout("input.temp");

  while(!filein.eof()){
    string temp;
    filein >> temp;
    if(temp==current) fileout << new_val << endl;
    else fileout << temp << endl;
  }
  rename("input.temp", "input.dat");
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BLOCKING METHOD ~~~~~~~~~~~~~~~~~~~~

void blocking_on_MD(int M, int N, string filename){

  double r;
  ofstream out;
  int L = M/N;
  out.open(filename, ios::out | ios::trunc);
  vector<double> sum(n_props,0);        // n_props index for n_props-dim measures
  vector<double> sum2(n_props,0);

  cout << "Starting simulation with blocking. " << endl;
  for(unsigned int i=0; i<N; i++){
    if((i)%5==0) cout << "Running block " << i << " of " << N << endl;
    vector<double> meas(n_props,0);
    for(int k=0; k<L; k++){
      Move();

      // istant print only of the 100-th measure
      if(k%100==0)Measure(true);
      else Measure(false);

      meas.at(iv)+=stima_pot;
      meas.at(ik)+=stima_kin;
      meas.at(ie)+=stima_etot;
      meas.at(it)+=stima_temp;
      // here i need Nbin individual operations of blocking meth
      // because each bin is a specific and different measure
      for(int ibin=4; ibin<n_props; ibin++){
          r = (ibin-4)*bin_size;
          stima_gofr[ibin-4] /= rho * npart * (((4 * pi )/ 3)* (pow(r + bin_size, 3) - pow(r, 3)));
          meas.at(ibin) += stima_gofr[ibin-4];
      }

    }

    for(int j=0; j<n_props; j++){
      sum.at(j) += meas.at(j)/(L/1.);
      sum2.at(j)+= (meas.at(j)/(L/1.))*(meas.at(j)/(L/1.));
      // need setprecision for showing the Etot fluctuations
      out << setprecision(8)  << sum.at(j)/(i+1) << "   " << setprecision(8) << error(sum.at(j)/(i+1), sum2.at(j)/(i+1), i) << "   ";
    }
    out << endl;
  }
  out.close();
  return;
}


double error(double val, double val2, unsigned int k){
  if(k==0) return 0;
  else return sqrt((val2-val*val)/k);
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOVE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RESCALING VELOCITIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ at first, accumulate a certain number of mean_v2 ~~~~~~~~~~~~~~~
// AGAIN, SAME COMMENTS DONE FOR LAB_04

double eval_mean_v2(){
  double mean_v2 = 0;

  for(int i=0; i<npart; i++){
    vx[i] = Pbc(x[i] - xold[i])/(delta);
    vy[i] = Pbc(y[i] - yold[i])/(delta);
    vz[i] = Pbc(z[i] - zold[i])/(delta);
    mean_v2+=vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }

  mean_v2 /= (double)npart;
  return mean_v2;
}

// ~~~~~~~~~~~~~ then use them for calculating a more stable value of sumv2 ~~~~~~~~~~~~~~~~~~~~~~~~~~

void rescale_velocities(vector<double> history_kin){

  double fs       = 0.0;
  double sumv2    = 0.0;

  for(unsigned int i=0; i<history_kin.size(); i++) sumv2 += history_kin.at(i);

  sumv2 /= history_kin.size();
  cout << "Media delle velocità medie calcolate: " << sumv2 << endl;

  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor

  cout << "Scale factor is " << fs << endl;
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FORCE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MEASURE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Measure(bool print_istant){ //Properties measurement
  int bin;
  long double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;


  if(print_istant==true){

    Epot.open("output_epot.dat",ios::app);
    Ekin.open("output_ekin.dat",ios::app);
    Temp.open("output_temp.dat",ios::app);
    Etot.open("output_etot.dat",ios::app);
  }

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

    bin = int(dr/bin_size);
    if(bin < nbins) stima_gofr[bin] += 2;     // into the hist

    if(dr < rcut){
      vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  stima_pot  = (v/(double)npart); //Potential energy per particle
  stima_kin  = (t/(double)npart); //Kinetic energy per particle
  stima_temp = ((2.0 / 3.0) * t/(double)npart); //Temperature
  stima_etot = ((t+v)/(long double)npart); //Total energy per particle

  if(print_istant==true){
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
  }
  return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ WRITE FINAL CONFIGURATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ AND FINAL t-dt CONFIGURATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteOld;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  WriteOld.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box    << "   " <<  y[i]/box    << "   " << z[i]/box << endl;
    WriteOld  << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  WriteOld.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAVE CONFIGUATION IN XYZ FORMAT ~~~~~~~~~~~~~~~~~~~~~~~~~

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
