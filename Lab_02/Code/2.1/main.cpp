#include <iostream>
#include "integrate.h"
#include "functions.h"
#include "statistic.h"
#include "random.h"
#include "measure.h"

int main(int argc, char* argv[]){

  if(argc < 4){
    cout << "You have to enter the number of throws M, the number of blocks N and the number of extraction for evaluating the integral!" << endl;
    return -1;
  }

  Random rnd;                                // classe Random
  Statistic mystat;                          // classe Statistic
  BaseFunction *f = new MyCos();             // concretizzo una FunzioneBase
                                             // nella forma della mia integranda
  // definisco una misura concreta di tipo Integrate
  Measure *I = new Integrate(0,1,atoi(argv[3]), f, "uniform sampling", rnd);
  // faccio blocking utilizzando I come misura all'interno
  mystat.blocking(atoi(argv[1]),atoi(argv[2]),I, "../../Results/uniform_sampling.dat");

  // concretizzo la FunzioneBase in un'altra funzione, utile per importance sampling
  f = new MyGeneral();
  // cambio il tipo di misura da usare in blocking
  I = new Integrate(0,1,atoi(argv[3]), f, "importance sampling", rnd);
  // faccio blocking utilizzando la nuova I come misura
  mystat.blocking(atoi(argv[1]),atoi(argv[2]),I, "../../Results/importance_sampling.dat");

  return 0;
}
