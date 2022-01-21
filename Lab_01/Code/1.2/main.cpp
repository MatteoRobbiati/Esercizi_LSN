#include <iostream>
#include "statistic.h"
#include "random.h"
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  int M = 1e4;
  Statistic mystat;
  vector<double> par_uni, par_exp, par_cau;
  par_uni.push_back(0);
  par_uni.push_back(1);
  par_exp.push_back(1);

  mystat.TLC(M, "uniform"    , par_uni, "../../Results/TLC_uniform.dat");
  mystat.TLC(M, "exponential", par_exp, "../../Results/TLC_exponential.dat");
  mystat.TLC(M, "cauchy"     , par_cau, "../../Results/TLC_cauchy.dat");

  return 0;
}
