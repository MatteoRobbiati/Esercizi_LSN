#include <iostream>
#include "econophysics.h"
#include "statistic.h"
#include "random.h"
#include "measure.h"

int main(int argc, char* argv[]){


  int M = 1e4;
  int N = 100;
  double asset = 100;
  double strike = 100;
  double period = 1;
  double rate = 0.1;
  double volatility = 0.25;
  string method = "direct";

  Random rnd;
  rnd.Init();
  Statistic mystat;
  Measure *first_option = new EconoPhysics(asset, period, rate, strike, volatility, method, rnd);
  mystat.blocking(M, N, first_option, "../../Results/direct_call.dat");

  return 0;
}
