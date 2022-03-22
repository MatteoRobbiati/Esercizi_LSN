#include <iostream>
#include "randomwalk.h"
#include "statistic.h"
#include "random.h"
#include "measure.h"

int main(int argc, char* argv[]){

  if(argc < 3){
    cout << "You have to enter the lattice constant and the number of steps of the walker!" << endl;
    return -1;
  }

  Random rnd;
  rnd.Init();

  int M = 1e5;
  double d = atoi(argv[1]);
  int N = atoi(argv[2]);
  int dimension = 3;


  RandomWalk disc_RW(d, 500, dimension, "discrete", rnd);
  RandomWalk cont_RW(d, 500, dimension, "continuum", rnd);

  disc_RW.run_walk("../../Results/discrete_RW.dat");
  cont_RW.run_walk("../../Results/continuum_RW.dat");

  Statistic mystat;
  Measure *rw_measure = new RandomWalk(d, N, dimension, "discrete", rnd);
  mystat.blocking(M, 100, rw_measure, "../../Results/block_on_steps_disc.dat");

  rw_measure = new RandomWalk(d,N,dimension, "continuum", rnd);
  mystat.blocking(M, 100, rw_measure, "../../Results/block_on_steps_cont.dat");


  return 0;
}
