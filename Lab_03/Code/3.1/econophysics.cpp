#include "econophysics.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;

EconoPhysics::EconoPhysics(double asset, double period, double rate, double strike, double volatility, string method, Random rnd){

  _rnd = rnd;
  _rnd.Init();
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
    params.push_back(0.0);
    params.push_back(t);
    double W = _rnd.generate_by_inversion("gaussian", params);
    return _asset*exp((_r-0.5*_volatility*_volatility)*t+_volatility*W);
  }
  else return 0;
}

vector<double> EconoPhysics::european_options(double t){
  vector<double> opt;
  if(_method=="direct"){
    double price = price_at(t);
    opt.push_back(exp(-_r*t)*max(0., price-_K));
    opt.push_back(exp(-_r*t)*max(0., _K-price));
  }
  if(_method=="step by step"){
    double price = _asset;
    double dt = t/100.;
    for(int i=0; i<100; i++){
      price = price*exp((_r-0.5*_volatility*_volatility)*dt+_volatility*_rnd.Gauss(0.,1.)*sqrt(dt));
    }
    opt.push_back(exp(-_r*t)*max(0., price-_K));
    opt.push_back(exp(-_r*t)*max(0., _K-price));
  }
  return opt;
}

double EconoPhysics::N(double d){
  return 0.5*(1.+erf(d/sqrt(2.)));
}

void EconoPhysics::save_Black_Scholes_solution(string filename){
  ofstream out;
  out.open(filename);

  double d1 = 1./(_volatility * sqrt(_T)) * (log(_asset / _K) + (_r + (_volatility*_volatility) / 2.) * _T);
  double d2 = (d1-_volatility*sqrt(_T));
  out << _asset*N(d1)-_K*exp(-_r*_T)*N(d2) << "   ";
  out << _asset*(N(d1)-1.)-_K*exp(-_r*_T)*(N(d2)-1.) << endl;
  out.close();

  return;
}


// ~~~~~~~~~~~~~~ concrete implementation of measure needed in blocking ~~~~~~~~~~~~~~~~~~~~~~~~~~

vector<double> EconoPhysics::get_measure(){     
    vector<double> opt = european_options(1.0);
    return opt;
}


int EconoPhysics::get_dimension(){
  int dim = 2;
  return dim;
}
