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

  Random rnd;
  Statistic mystat;
  BaseFunction *f = new MyCos();
  Measure *I = new Integrate(0,1,atoi(argv[3]), f, "uniform sampling", rnd);

  mystat.blocking(atoi(argv[1]),atoi(argv[2]),I, "../../Results/uniform_sampling.dat");

  f = new MyGeneral();
  I = new Integrate(0,1,atoi(argv[3]), f, "importance sampling", rnd);

  mystat.blocking(atoi(argv[1]),atoi(argv[2]),I, "../../Results/importance_sampling.dat");

  return 0;
}
