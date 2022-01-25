#include "integrate.h"
#include <string>
#include <iostream>

using namespace std;

Integrate::Integrate(double a, double b, int n, BaseFunction *integrating, string integration, Random rnd){
  _rnd = rnd;
  _rnd.Init();

  _g = integrating;
  _n = n;
  _type = integration;
  _a = a;
  _b = b;

  _dim = 1;

  if(_b>_a) _sign = 1;
  else _sign = -1;

  if(_type!="uniform sampling" && _type!="importance sampling"){
    cout << "You must chose one of these type of integration: uniform sampling, importance sampling." << endl;
  }
}

Integrate::~Integrate(){
  _rnd.SaveSeed();
}

vector<double> Integrate::get_measure(){
    vector<double> meas;
    double sum = 0;

    if(_type=="uniform sampling"){
      for(int i=0; i<_n; i++) sum+=_g->eval(_rnd.Rannyu());
    }
    if(_type=="importance sampling"){
      for(int i=0; i<_n; i++) sum+=_g->eval(1.+sqrt(1.-_rnd.Rannyu()));
    }
    meas.push_back(_sign*sum/_n);
    return meas;
}

int Integrate::get_dimension(){
  return _dim;
}
