#include "thismeasure.h"
#include "measure.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

ThisMeasure :: ThisMeasure(Random rnd){
  _dim = 2;                                       // valore e scarto quadratico in questo caso
  _rnd = rnd;
  _rnd.Init();
};

ThisMeasure ::~ThisMeasure(){
  _rnd.SaveSeed();
};

int ThisMeasure :: get_dimension(){
  return _dim;
}

vector<double> ThisMeasure :: get_measure(){
  vector<double> meas;                            // qui implemento la misura
  double r = _rnd.Rannyu();
  meas.push_back(r);
  meas.push_back((r-0.5)*(r-0.5));
  return meas;
}
