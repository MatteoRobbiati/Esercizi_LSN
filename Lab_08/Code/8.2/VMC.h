/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __VMC__
#define __VMC__

//Random numbers
#include "random.h"
#include <string>
#include <iostream>

using namespace std;

int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,ie;
double nbins;
double stepsize;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_e;
double block_e;                     // here i save the last estimation in a run
double err_e;

// thermodynamical state
double beta,temp, mu, sigma, delta;

// variable 1d for sampling

double x;

// simulation
int nstep, nblk, metro, L, print_hist, explore;
string restart;

//functions
void Read_input();
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move();
void Explore();
void Measure(bool);
void Single_run(double mu, double sigma);
double Error(double,double,int);
void set_restart(string, string);
void Equilibrate_system(int);

void print_x(string);
double psi(double x, double mu, double sigma);
double gauss(double x, double mu, double sigma);
double eval_V();
double apply_H(double mu, double sigma);
double p(double x);
void show_exploration(string);
void simulated_annealing(double betai, double betaf, int Nstep);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
