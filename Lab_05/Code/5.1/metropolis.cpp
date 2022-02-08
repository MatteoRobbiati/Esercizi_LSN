#include "metropolis.h"

use namespace std;

Metropolis::Metropolis(int N, pdf *mypdf, Position *start, Random rnd){
  _steps = N;
  _accepted = 0;

  _x = start;
  _xold->to_zero();

  p_function = mypdf;

  _rnd = rnd;
  rnd.Init();
}

Metropolis::~Metropolis(){
  _rnd.SaveSeed();
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double Metropolis::T_uni(){
  
}
