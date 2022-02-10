#include <iostream>
#include "metropolis.h"
#include "pdf.h"
#include "position.h"
#include "statistic.h"
#include <string>
#include "random.h"
#include "measure.h"

using namespace std;

int main(int argc, char* argv[]){

  // some useful variables

  int M = 1e6;
  double step_100 = 1.235;
  double step_210 = 2.955;

  Random* rnd = new Random();
  rnd->Init();
  Statistic *mystat = new Statistic();

  // some useful positions in the 3d space; Position is a specific class

  Position *start  = new Position(1.,1.,1., rnd);
  vector<double> coord  = start->get_coordinates();           // saving start position
  vector<double> origin(3,0.);

  // Simulating the orbitals: pdf is a virtual class whose implementation are hydro_100 and hydro_210

  cout << "Running Metropolis algo for simulating |100> and |210> states of hydrogen. " << endl;

  pdf *hydrogen = new hydro_100();
  Metropolis *Metro_100 = new Metropolis(M, hydrogen, start, rnd, step_100, "uniform");
  Metro_100->run("../../Results/H100.dat");
  cout << "rate hydrogen |100>: " << Metro_100->rate_of_acceptance() << endl;

  start->set_coordinates(coord);                              // starting from the old start

  hydrogen = new hydro_210();
  Metropolis *Metro_210 = new Metropolis(M, hydrogen, start, rnd, step_210, "uniform");
  Metro_210->run("../../Results/H210.dat");
  cout << "rate hydrogen |210>: " << Metro_210->rate_of_acceptance() << endl;

  // Evaluation of the expected values

  cout << "Evaluating expected values for |Psi_100|^2 and |Psi_210|^2 with a blocking method." << endl;
  cout << endl;
  start->set_coordinates(coord);

  hydrogen = new hydro_100();
  Measure *measure_100 = new Metropolis(M, hydrogen, start, rnd, step_100, "gaussian");
  mystat->blocking(M, 100, measure_100, "../../Results/gau_blockin_on_100.dat");

  cout << endl;
  start->set_coordinates(coord);

  hydrogen = new hydro_210();
  Measure *measure_210 = new Metropolis(M, hydrogen, start, rnd, step_210, "gaussian");
  mystat->blocking(M, 100, measure_210, "../../Results/gau_blockin_on_210.dat");



  return 0;
}
