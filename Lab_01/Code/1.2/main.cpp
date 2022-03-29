#include <iostream>
#include "statistic.h"
#include "random.h"
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  int M = 1e4;
  Statistic mystat;
  vector<double> par_uni, par_exp, par_cau;
  // ho scritto il metodo per generare variabili con Box-Muller in modo generico
  // perciò serve fornire un vettore di parametri che definiscano la distribuzione
  // vector generico perché ogni distribuzione ha un num variabile di parametri
  par_uni.push_back(0);
  par_uni.push_back(1);
  par_exp.push_back(1);

  mystat.TLC(M, "uniform"    , par_uni, "../../Results/TLC_uniform.dat");
  mystat.TLC(M, "exponential", par_exp, "../../Results/TLC_exponential.dat");
  mystat.TLC(M, "cauchy"     , par_cau, "../../Results/TLC_cauchy.dat");

  return 0;
}
