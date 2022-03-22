#include <iostream>
#include "statistic.h"
#include "random.h"
#include "measure.h"
#include "thismeasure.h"
#include <iomanip>
#include <stdio.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  int M = 1e4;
  int N = 1e2;
  ofstream out;
  Random rnd;
  Statistic mystat;
  Measure *check_uniform = new ThisMeasure(rnd);


  // PARTE 1
  // il metodo blocking è implementato in statistic.cpp.
  // è composto da block_step, in cui genero le misure A[] (medie su L realizzazioni)
  // e da blocking, cioè la media iterativa sulle misure A
  mystat.blocking(M, N, check_uniform, "../../Results/results.dat");

  // PARTE 2
  // N = 100 volte
  // calcolo del chi quadro: genero M valori e divido lo spazio campione in N sezioni
  // calcolo chiquad usando i campionamenti di ogni sezione messi a confronto con il valore atteso

  out.open("../../Results/chiquad.dat", ios::out | ios::trunc);
  for(int j=0; j<100; j++){
    vector<double> counter(N,0);
    for(int i=0; i<M; i++){
      int rand = int(mystat.uniform_sampling(0,1)*100);
      counter.at(rand)++;
    }
    out << setprecision(6) << mystat.chiquad(counter, M/N, true) << endl;        // stampo così da plottare
  }                                                           // l'istogramma

  out.close();

  return 0;
}
