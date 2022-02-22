#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "random.h"
#include "TSP.h"

using namespace std;

int main(){

  int N = 32;
  int M = 1000;
  string circle_circuit = "../../Results/circle_circuit.dat";
  string big_circle_cir = "../../Results/bigcircle.dat";
  string square_circuit = "../../Results/square_circuit.dat";
  Random rnd;
  rnd.Init();
  rnd.generate_cities_on_circle(1, N, circle_circuit);
  rnd.generate_cities_on_circle(1, 100, big_circle_cir);
  rnd.generate_cities_in_square(1, N, square_circuit);

  Salesman *Jim_on_circle = new Salesman(N, M, circle_circuit, rnd);
  Salesman *Jim_in_square = new Salesman(N, M, square_circuit, rnd);
  Salesman *Jim_on_CIRCLE = new Salesman(100, M, big_circle_cir, rnd);

  cout << "Jim is now working on the cities displaced on the circle! " << endl;

  Jim_on_circle->show_best_chromo();
  Jim_on_circle->save_chromo(0,"../../Results/circle_initial.dat", true);
  Jim_on_circle->run(1000);
  Jim_on_circle->save_chromo(0,"../../Results/circle_final.dat", true);
  Jim_on_circle->show_best_chromo();

  cout << endl << "Jim is now working into the square!" << endl;

  Jim_in_square->show_best_chromo();
  Jim_in_square->save_chromo(0,"../../Results/square_initial.dat", true);
  Jim_in_square->run(1000);
  Jim_in_square->save_chromo(0,"../../Results/square_final.dat", true);
  Jim_in_square->show_best_chromo();

  cout << endl << "Jim is now working on the cities displaced on the BIG circle! " << endl;

  Jim_on_CIRCLE->show_best_chromo();
  Jim_on_CIRCLE->save_chromo(0,"../../Results/CIRCLE_initial.dat", true);
  Jim_on_CIRCLE->run(10000);
  Jim_on_CIRCLE->save_chromo(0,"../../Results/CIRCLE_final.dat", true);
  Jim_on_CIRCLE->show_best_chromo();

  return 0;
}
