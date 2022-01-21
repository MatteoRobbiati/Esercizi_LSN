#include <iostream>
#include "statistic.h"
#include "random.h"
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  int M = 1e4;
  int N = 1e2;
  vector<double> counter(N,0);
  Statistic mystat;

  // PARTE 1
  // il metodo blocking è implementato in statistic.cpp.
  // è composto da block_step, in cui genero le misure A[] (medie su L realizzazioni)
  // e da blocking, cioè la media iterativa sulle misure A
  mystat.blocking(M, N, "../../Results/results.dat");

  // PARTE 2
  // calcolo del chi quadro: genero M valori e divido lo spazio campione in N sezioni
  // calcolo chiquad usando i campionamenti di ogni sezione messi a confronto con il valore atteso
  for(int i=0; i<M; i++){
    int rand = int(mystat.uniform_sampling(0,1)*100);
    counter.at(rand)++;
  }

  cout << "Il chi quadro calcolato è: " << mystat.chiquad(counter, M/N, true) << endl;
  return 0;
}
