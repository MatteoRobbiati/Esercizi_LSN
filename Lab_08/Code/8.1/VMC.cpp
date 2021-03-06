/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include "VMC.h"
#include <string>

using namespace std;

int main(){

  Read_input();                // input
  Equilibrate_system(10000);   // equilibration
  Single_run(mu, sigma);       // a single run of blocking

  return 0;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  A SINGLE RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Single_run(double mu, double sigma){
  for(int iblk=1; iblk <= nblk; ++iblk){            //Simulation
    Reset(iblk);                                    //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move();
      Measure();
      Accumulate();                                 //Update block averages
    }
    Averages(iblk);                                 //Print results for current block
  }
  return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PREPARATION OF THE ALGO ~~~~~~~~~~~~~~~~~~~~~~~~~~

void Read_input(){

  n_props = 1;
  ie = 0;           // total energy index for my vecs

  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  ifstream in;
  in.open("input.file");

  in >> x;
  in >> stepsize;
  in >> print_hist;
  in >> mu;
  in >> sigma;
  in >> restart;
  in >> explore;
  in >> delta;
  in >> nblk;
  in >> nstep;

  // some generalities of the simulation
  if(explore==0){
    cout << "Code is running with a single evaluation of the energy. " << endl;
    cout << "It starts with the choice: (mu, sigma) = (" << mu << ", " << sigma << ")" << endl;
  }else{
    cout << "Code is running exploring the parameter space." << endl;
    cout << "It will move a walker with steps of modulus " << delta << " in the space.";
  }

  cout << "--------------------------------------------------------------------" << endl;
  cout << "Other generalities of the algo: " << endl;
  cout << "starting  x = " << x << endl;
  cout << "stepsize dx = " << stepsize << endl;
  cout << "restart opt = " << restart << endl;
  cout << "# of blocks = " << nblk << endl;
  cout << "# of steps  = " << nstep << endl;

  in.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SWAP RESTART VALUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void set_restart(string current, string new_val){
  string strReplace = current;
  string strNew = new_val;
  ifstream filein("input.file");
  ofstream fileout("input.temp");

  while(!filein.eof()){
    string temp;
    filein >> temp;
    if(temp==current) fileout << new_val << endl;
    else fileout << temp << endl;
  }
  rename("input.temp", "input.file");
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EQUILIBRATION FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Equilibrate_system(int equi_steps){
  cout << "Restart value is " << restart << endl;
  if(restart=="true"){
    cout << "Equilibration phase is running, walker will be moved " << equi_steps << " times" << endl;
    cout << "And those values will be ignored." << endl;   // Moving for equi_steps-time the walker
    for(int i=0; i<equi_steps; i++) Move();                // in the ending of EQUILIBRATION
    set_restart("true","false");                           // change restart value true->false
  }else{
    cout << "No equilibration phase needed " << endl;
  }
  return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOVE FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Move(){
  // metropolis uses squared modulus of the wave-function
  double x_new = x + rnd.Rannyu(-stepsize, stepsize);
  double alpha = min(1., p(x_new)/p(x));
  double ran   = rnd.Rannyu();

  if(ran < alpha){ x = x_new; accepted++;};
  attempted++;
  if(print_hist==1) print_x("sampling.dat");

  return;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MEASURE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Measure(){
  ofstream out;
  out.open("energia_one_shot.dat", ios::app);
  walker[ie] = apply_H(mu, sigma) + eval_V();
  out << walker[ie] << endl;
  out.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// gaussian
double gauss(double x, double mu, double sigma){
  return exp(-(x-mu)*(x-mu)/(2.*sigma*sigma));
}

// WF
double psi(double x, double mu, double sigma){
  return gauss(x, mu, sigma)+gauss(x, -mu, sigma);
}

// action of H on WF in x
double apply_H(double mu, double sigma){
  double kin = 1./(2*sigma*sigma)*(1.-1./psi(x, mu, sigma)*((x-mu)*(x-mu)/(sigma*sigma)*gauss(x,mu,sigma)+(x+mu)*(x+mu)/(sigma*sigma)*gauss(x,-mu,sigma)));
  return kin;
}

// potential in x
double eval_V(){
  return (x*x-2.5)*x*x;
}

// the p.d.f.
double p(double x){
  return psi(x, mu, sigma)*psi(x, mu, sigma);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SOME USEFUL METHODS FOR BLOCKING ~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Reset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Reset(int iblk){
   if(iblk == 1){
       for(int i=0; i<n_props; ++i){
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }
   for(int i=0; i<n_props; ++i) blk_av[i] = 0;
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Accumulate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Accumulate(void){
   for(int i=0; i<n_props; ++i) blk_av[i] = blk_av[i] + walker[i];
   blk_norm = blk_norm + 1.0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Averages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Averages(int iblk){

   ofstream Ene;
   const int wd=12;

    if(explore==0){
      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << accepted/double(attempted) << endl << endl;
    }

    Ene.open("ene.dat",ios::app);
    stima_e = blk_av[ie]/blk_norm; //Energy
    glob_av[ie]  += stima_e;
    glob_av2[ie] += stima_e*stima_e;
    err_e=Error(glob_av[ie],glob_av2[ie],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_e << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_e << endl;

    Ene.close();
    return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double Error(double sum, double sum2, int iblk){
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Print x value on a file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// for the histogram
void print_x(string filename){
  ofstream out;
  out.open(filename, ios::app);
  out << x << endl;
  out.close();
  return;
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
