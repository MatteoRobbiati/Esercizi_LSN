#include "randomwalk.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

// ~~~~~~~~~~~  constructor and destructor ~~~~~~~~~~~~~~~~~~~~

RandomWalk :: RandomWalk(double d, int N, int dim, string metric, Random rnd){
  _d = d;
  _N = N;
  _dim = dim;
  _metric = metric;

  for(int i=0; i<_dim; i++) _pos.push_back(0.);

  _rnd = rnd;
  _rnd.Init();
}

RandomWalk ::~RandomWalk(){
  _rnd.SaveSeed();
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~ discrete metric ~~~~~~~~~~~~~~~~~~~~~~~

void RandomWalk::make_discrete_step(){
  double rand = _rnd.Rannyu(-_dim,_dim);
  int versus = rand/abs(rand);
  int index = versus*int(rand);
  _pos.at(index)+=versus*_d;
  return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~ continuum metric ~~~~~~~~~~~~~~~~~~~~~~~

void RandomWalk::make_continuum_step(){
  double phi =   _rnd.Rannyu(0,2*M_PI);
  double theta = _rnd.Rannyu(0,M_PI);
  _pos.at(0) += _d*cos(phi)*sin(theta);
  _pos.at(1) += _d*sin(phi)*sin(theta);
  _pos.at(2) += _d*cos(theta);
  return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~ one random walk ~~~~~~~~~~~~~~~~~~~~~~~~~~~

void RandomWalk::run_walk(string filename){
  ofstream out;
  out.open(filename);
  if(_metric=="discrete"){
    for(int i=0; i<_N; i++){
      make_discrete_step();
      for(int j=0; j<_dim; j++) out << _pos.at(j) << "    ";
      out << endl;
    }
    return;
  }if(_metric=="continuum"){
    for(int i=0; i<_N; i++){
      make_continuum_step();
      for(int j=0; j<_dim; j++) out << _pos.at(j) << "    ";
      out << endl;
    }
    return;
  }
  else cout << "Please, choose one metric: discrete or continuum." << endl;
  out.close();
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~ some utilities ~~~~~~~~~~~~~~~

vector<double> RandomWalk::get_position(){
  return _pos;
}

void RandomWalk::walker_at_home(){
  for(int i=0; i<_dim; i++) _pos.at(i)=0;
}


double RandomWalk::get_squared_mod(){
  return _pos.at(0)*_pos.at(0)+_pos.at(1)*_pos.at(1)+_pos.at(2)*_pos.at(2);
}

// ~~~~~~~~~~~~~~~~~~~~~~~ concrete implementation of measure ~~~~~~~~

vector<double> RandomWalk::get_measure(){

  vector<double> step_distance;       // in each slot the sum of _N dist reached at that step
  walker_at_home();
  step_distance.push_back(get_squared_mod());

  if(_metric=="discrete"){
    for(int step=0; step<_N; step++){
      make_discrete_step();
      step_distance.push_back(get_squared_mod());
    }
  }if(_metric=="continuum"){
    for(int step=0; step<_N; step++){
      make_continuum_step();
      step_distance.push_back(get_squared_mod());
    }
  }

  return step_distance;              // la misura restituisce un vector
}

int RandomWalk::get_dimension(){
  return _N;                         // la misura contiene _N elementi
}
