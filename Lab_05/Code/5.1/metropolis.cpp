#include "metropolis.h"
#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

Metropolis::Metropolis(int N, pdf *mypdf, Position *start, Random* rnd, double stepsize, string method){
  _steps = N;
  _accepted = 0;
  _attempted=0;
  _stepsize = stepsize;
  _method = method;

  _x = start;

  p_function = mypdf;

  _rnd = rnd;
}

Metropolis::~Metropolis(){
  _rnd->SaveSeed();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EQUILIBRATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Metropolis::Equilibrate(int Nequi){
  ofstream out;
  out.open("../../Results/equilibration.dat");

  for(int i=0; i<Nequi; i++){
    try_step();
    out << _x->get_radius() << endl;
  }

  out.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRY STEP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Metropolis::try_step(){

  vector<double> old = _x->get_coordinates();
  double A =  p_function->eval(_x);

  if(_method=="uniform")   _x->uniform_step(_stepsize);
  if(_method=="gaussian")  _x->gaussian_step(_stepsize);

  double B =  p_function->eval(_x);
  double alpha = min(1., B/A);

  double rand = _rnd->Rannyu();
  if(rand<=alpha){
    _accepted++;
  }else _x->set_coordinates(old);
  _attempted++;
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ A RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Metropolis::run(string filename){
  ofstream out;
  out.open(filename);

  for(int i=0; i<_steps; i++){
    try_step();
    vector<double> coord = _x->get_coordinates();
    out << coord.at(0) << "   " << coord.at(1) << "   " << coord.at(2) << endl;
    //_x->show_position();
  }
  out.close();
  return;
}

double Metropolis::rate_of_acceptance(){
  return double(_accepted)/double(_attempted);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~ concrete implementation of measure ~~~~~~~~~~~~~~~~~~~~~

vector<double> Metropolis::get_measure(){
  try_step();
  vector<double> meas;
  meas.push_back(_x->get_radius());
  return meas;
}

int Metropolis::get_dimension(){
  return 1;
}
