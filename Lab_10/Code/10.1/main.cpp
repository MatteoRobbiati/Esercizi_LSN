#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "random.h"
#include "TSP.h"

using namespace std;

int main(){

  int N = 32;
  string circle_circuit = "../../Results/circle_circuit.dat";

  Random rnd;
  rnd.Init();
  rnd.generate_cities_on_circle(1, N, circle_circuit);

  Salesman *Jim_on_circle = new Salesman(N, circle_circuit, rnd);
  cout << "Jim is now working on the cities displaced on the circle! " << endl;
  Jim_on_circle->show_chromo();
  Jim_on_circle->save_chromo("../../Results/circle_ini.dat", true);
  Jim_on_circle->simulated_annealing(1., 1000, 1000, 200);
  Jim_on_circle->show_chromo();
  Jim_on_circle->save_chromo("../../Results/circle_fin.dat", true);


  return 0;
}
