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
#include "Monte_Carlo_ISING_1D.h"
#include <string>

using namespace std;

int main(){

  Input("T=inf", false);
  Equilibrate_system(1000);

  for(int iT = 0; iT<=((2.-0.5)/stepT)*varyTemp; iT++){    // only 1 loop if varyTemp is 0
    if(iT==0) cout << "Starting simulation at T=" << temp << endl << "----------------------------------" << endl;
    if(iT%1==0) cout << "Running evaluation for T=" << temp << endl;
    beta = 1./temp;
    for(int iblk=1; iblk <= nblk; ++iblk){            //Simulation
      Reset(iblk);                                    //Reset block averages
      print_flag = 0;
      for(int istep=1; istep <= nstep; ++istep){
        Move(metro);
        Measure();
        Accumulate();                                 //Update block averages
        print_flag++;
      }
      Averages(iblk);                                 //Print results for current block
    }
    temp += stepT;
  }
  ConfFinal();
  return 0;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INPUT FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Input(string initial_condition, bool print){
  ifstream ReadInput;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.file");

  ReadInput >> temp;
  ReadInput >> nspin;
  ReadInput >> J;
  ReadInput >> h;
  ReadInput >> metro;            // if=1 Metropolis else Gibbs
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> restart;
  ReadInput >> varyTemp;
  ReadInput >> stepT;
  L = nstep/nblk;           //steps per block

// if T is fixed i define beta as 1/T, if T is wanted to vary i initialize T at 0.5, the extreme of the range [0.5,2]
  if (varyTemp == 0){
    beta = 1./temp;
  }else{
      if(varyTemp == 1){
        cout << "Temperature will vary from 0.5 to 2.0 during the simulation" << endl;
        temp = 0.5;
        beta = 1./temp;
      } else {
          cerr << " *** Wrong varyTemp value: it should be 0 or 1 ***" << endl;
          exit(1);
        }
  }

  if(print==true){
    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;
    cout << "Temperature = " << temp << endl;
    cout << "Number of spins = " << nspin << endl;
    cout << "Exchange interaction = " << J << endl;
    cout << "External field = " << h << endl << endl;
    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
  }
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  n_props = 4; //Number of observables

// read initial config
// from config.0 if restart is true, from config.final if it is false

  if(restart=="true"){
    prepare_initial_config(initial_condition, "config.0");
    ifstream in;
    in.open("config.0");
    for(int i=0; i<nspin; i++) in >> s[i];
    in.close();
  }else{
    ifstream in;
    in.open("config.final");
    for(int i=0; i<nspin; i++) in >> s[i];
    in.close();
  }

//Evaluate energy etc. of the initial configuration
  Measure();
//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HOW TO PREPARE CONFIG.0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// all +1 if T=0 and randomly +/- 1 if T=inf
void prepare_initial_config(string init, string filename){
  ofstream out;
  out.open(filename);

  if(init == "T=0"){
    cout << "Initial condition is T=0" << endl;
    for(int i=0; i<nspin; ++i){
       s[i] = 1;
       out << s[i] << endl;
    }
  }
  else{
    for(int i=0; i<nspin; ++i){
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
      out << s[i] << endl;
    }
  }
  out.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SWAP RESTART VALUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EQUILIBRATION FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Equilibrate_system(int equi_steps){
  cout << "Restart value is " << restart << endl;
  cout << "Equilibration phase is running, walker will be moved " << equi_steps << " times and final configuration ";
  cout << "of the spins will be saved in file config.final" << endl;
  for(int i=0; i<equi_steps; i++) Move(metro);           // Moving for equi_steps-time the walker
  set_restart("true","false");                           // in the ending of EQUILIBRATION
                                                         // change restart value true->false
  ConfFinal();
  return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOVE FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Move(int metro)
{
  int o;
  double deltaE, r, alpha, p;
  //double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;

  for(int i=0; i<nspin; ++i){
    o = (int)(rnd.Rannyu()*nspin);
    if(metro==1){                                               //Metropolis
      deltaE = -2 * Boltzmann(s[o], o);
      alpha = min(1., exp(-beta*deltaE));
      attempted += 1;
      r = rnd.Rannyu();
      if (r <= alpha) {
        s[o] *= -1;
        accepted += 1;
      }
    }                                                          // Gibbs
    else{
      deltaE = -2 * Boltzmann(s[o], o) * s[o]; // s[o] * s[o] = 1
      p = 1./(1 + exp(-beta * deltaE)); // probability to have p( s_o = 1 | {s_j : j=/=o})
      attempted += 1;
      accepted += 1;
      // Decide what value should assume spin s_o
      r = rnd.Rannyu();
      if (r <= p) {s[o] = +1;}
      else  	    {s[o] = -1;}
    }
  }
  return;
}

double Boltzmann(int sm, int ip){
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;     // n.n. of o-th spin
  return ene;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MEASURE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Measure(){
//  int bin;
  long double u = 0.0, m = 0.0;
//cycle over spins
  for (int i=0; i<nspin; ++i){
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u * u;
  walker[im] = m;
  walker[ix] = beta * m * m;

  if(varyTemp==0){
    ofstream out;
    out.open("istant_u.dat", ios::app);
    out << u/(double)nspin << endl;
  }
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

   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;

    if(varyTemp==0){
      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }

    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    if(varyTemp==1){
      //cout << attempted << endl;
      if(iblk==20) Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    }else{
      Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    }
    Ene.close();

    Heat.open("output.heat.0", ios::app);
    stima_c = beta*beta*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double)nspin; // Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c = Error(glob_av[ic],glob_av2[ic],iblk);
    if(varyTemp==1){
      if(iblk==20) Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    }else{
      Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    }
    Heat.close();

    Mag.open("output.mag.0", ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; // Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m = Error(glob_av[im],glob_av2[im],iblk);
    if(varyTemp==1){
      if(iblk==20) Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    }else{
      Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    }
    Mag.close();

    Chi.open("output.chi.0", ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin; // Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x = Error(glob_av[ix],glob_av2[ix],iblk);
    if(varyTemp==1){
      if(iblk==20) Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    }else{
      Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    }
    Chi.close();
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PRINT FINAL CONFIG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void ConfFinal(void){
  ofstream WriteConf;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i) WriteConf << s[i] << endl;
  WriteConf.close();
  rnd.SaveSeed();
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PBC AND ERROR (FOR BLOCKING) ~~~~~~~~~~~~~~~~~~~~~~~~~

int Pbc(int i){  //Algorithm for periodic boundary conditions
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk){
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
