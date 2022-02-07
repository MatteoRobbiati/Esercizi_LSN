#include <iostream>
#include "econophysics.h"
#include "statistic.h"
#include "random.h"
#include "measure.h"


int main(int argc, char* argv[]){


  int M = 1e5;
  int N = 100;
  double asset = 100.;
  double strike = 100.;
  double period = 1.;
  double rate = 0.1;
  double volatility = 0.25;
  string method [] = {"direct", "step by step"};

  Random rnd;
  Statistic mystat;
  EconoPhysics MyEco(asset, period, rate, strike, volatility, method[0], rnd);

  //~~~~~~~~~~~~~~~~ analytic's solutions ~~~~~~~~~~~~~~~~~~~~~~~~~~

  MyEco.save_Black_Scholes_solution("../../Results/solutions.dat");

  // ~~~~~~~~~~~~~~ simulations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Measure *options = new EconoPhysics(asset, period, rate, strike, volatility, method[0], rnd);
  mystat.blocking(M, N, options, "../../Results/direct_opt.dat");

  options = new EconoPhysics(asset, period, rate, strike, volatility, method[1], rnd);
  mystat.blocking(M, N, options, "../../Results/bystep_opt.dat");

  return 0;
}
