#include <iostream>
#include "metropolis.h"
#include "pdf.h"
#include "position.h"
#include <string>
#include "random.h"
#include "measure.h"

using namespace std;

int main(int argc, char* argv[]){


  int M = 1e5;
  double step_100 = 1.235;
  double step_210 = 4.95;
  Random* rnd = new Random();
  rnd->Init();
  Position *start = new Position(100.,100.,100., rnd);
  vector<double> coord = start->get_coordinates();

  pdf *hydrogen = new hydro_100();
  Metropolis *MyMetro = new Metropolis(M, hydrogen, start, rnd, step_100, "uniform");
  MyMetro->run("../../Results/H100.dat");
  cout << "rate hydrogen |100>: " << MyMetro->rate_of_acceptance() << endl;

  start->set_coordinates(coord);

  hydrogen = new hydro_210();
  Metropolis *Metro_210 = new Metropolis(M, hydrogen, start, rnd, step_210, "uniform");
  Metro_210->run("../../Results/H210.dat");
  cout << "rate hydrogen |210>: " << Metro_210->rate_of_acceptance() << endl;

  return 0;
}
