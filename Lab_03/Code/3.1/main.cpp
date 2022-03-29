#include <iostream>
#include "econophysics.h"
#include "statistic.h"
#include "random.h"
#include "measure.h"


int main(int argc, char* argv[]){

  // blockParameters
  int M = 1e5;
  int N = 100;
  // EconoParameters
  double asset = 100.;
  double strike = 100.;
  double period = 1.;
  double rate = 0.1;
  double volatility = 0.25;
  // due metodi di simulazione
  string method [] = {"direct", "step by step"};

  Random rnd;
  Statistic mystat;
  // la classe EconoPhysics
  EconoPhysics MyEco(asset, period, rate, strike, volatility, method[0], rnd);

  //~~~~~~~~~~~~~~~~ analytic's solutions ~~~~~~~~~~~~~~~~~~~~~~~~~~

  // che uso per le soluzioni di B&S
  MyEco.save_Black_Scholes_solution("../../Results/solutions.dat");

  // ~~~~~~~~~~~~~~ simulations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // come misura concreta una EcoPhy con metodo diretto
  Measure *options = new EconoPhysics(asset, period, rate, strike, volatility, method[0], rnd);
  mystat.blocking(M, N, options, "../../Results/direct_opt.dat");
  // e una con metodo step by step
  options = new EconoPhysics(asset, period, rate, strike, volatility, method[1], rnd);
  mystat.blocking(M, N, options, "../../Results/bystep_opt.dat");

  return 0;
}
