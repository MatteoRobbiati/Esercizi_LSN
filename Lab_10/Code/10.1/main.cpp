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
  string square_circuit = "../../Results/square_circuit.dat";

  Random rnd;
  rnd.Init();
  //rnd.generate_cities_on_circle(1, N, circle_circuit);
  rnd.generate_cities_in_square(1, N, square_circuit);

  /*Salesman *Jim_on_circle = new Salesman(N, circle_circuit, rnd);
  cout << "Jim is now working on the cities displaced on the circle! " << endl;
  Jim_on_circle->show_chromo();
  Jim_on_circle->save_chromo("../../Results/circle_ini.dat", true);
  Jim_on_circle->simulated_annealing(1., 500, 500, 100, "../../Results/anneal_circle.dat");
  Jim_on_circle->show_chromo();
  Jim_on_circle->save_chromo("../../Results/circle_fin.dat", true);*/

  Salesman *Jim_in_square = new Salesman(N, square_circuit, rnd);
  cout << "Jim is now working on the cities displaced in the square! " << endl;
  Jim_in_square->show_chromo();
  Jim_in_square->save_chromo("../../Results/square_ini.dat", true);
  Jim_in_square->simulated_annealing(1., 500, 500, 100, "../../Results/anneal_square.dat");
  Jim_in_square->show_chromo();
  Jim_in_square->save_chromo("../../Results/square_fin.dat", true);


  return 0;
}
