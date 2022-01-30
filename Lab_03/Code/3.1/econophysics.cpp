#include "econophysics.h"
#include <cmath>
#include <iostream>

using namespace std;

EconoPhysics::EconoPhysics(double asset, double period, double rate, double strike, double volatility, string method, Random rnd){

  _rnd = rnd;
  _asset = asset;
  _T = period;
  _r = rate;
  _K = strike;
  _volatility = volatility;
  _method = method;
}

EconoPhysics::~EconoPhysics(){
  _rnd.SaveSeed();
}

// ~~~~~~~~~~~~~~~~~~~~~ Eco methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double EconoPhysics::price_at(double t){
  if(_method=="direct"){
    vector <double> params;
    params.push_back(0);
    params.push_back(t);
    double W = _rnd.generate_by_inversion("gaussian", params);
    return _asset*exp((_r-0.5*_volatility*_volatility)*t+_volatility*W);
  }
  else return 0;
}

double EconoPhysics::d1_at(double t){
  return (1./sqrt(_T-t))*(log(price_at(t)/_K)+(_r+0.5*_volatility*_volatility*(_T-t)));
}

double EconoPhysics::d2_at(double t){
  return d1_at(t)-_volatility*sqrt(_T-t);
}

double EconoPhysics::N(double d){
  return 0.5*(1+erf(d/sqrt(2)));
}

double EconoPhysics::call_option(double t){
  double d1 = d1_at(t);
  double d2 = d1-_volatility*sqrt(_T-t);
  cout << d1 << endl;
  return price_at(t)*N(d1)-_K*exp(-_r*(_T-t))*N(d2);
}




// ~~~~~~~~~~~~~~ concrete implementation of measure needed in blocking ~~~~~~~~~~~~~~~~~~~~~~~~~~

vector<double> EconoPhysics::get_measure(){
  vector<double> measure;

  if(_method=="direct"){
    measure.push_back(call_option(0.0));
    return measure;
  }

  else return measure;
}

int EconoPhysics::get_dimension(){
  int dim = 1;
  return dim;
}
