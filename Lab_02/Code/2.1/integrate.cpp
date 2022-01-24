#include "integrate.h"

using namespace std;

Integrate::Integrate(double a, double b, int n, BaseFunction *integrating, Random rnd){
  _rnd = rnd;
  _rnd.Init();

  _g = integrating;
  _a = a;
  _b = b;
  _n = n;

  _dim = 1;

  if(_b>_a) _sign = 1;
  else _sign = -1;
}

Integrate::~Integrate(){
  _rnd.SaveSeed();
}

vector<double> Integrate::get_measure(){
  vector<double> meas;
  double sum = 0;
  for(int i=0; i<_n; i++) sum+=_g->eval(_rnd.Rannyu());
  meas.push_back(_sign*sum/_n);
  return meas;
}

int Integrate::get_dimension(){
  return _dim;
}
