#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "random.h"
#include "TSP.h"

using namespace std;

int main(){

  ofstream out;
  out.open("../../Results/hist.dat");

  int N = 32;
  int M = 1000;
  string circuit = "../../Results/circle_circuit.dat";
  Random rnd;
  rnd.Init();
  rnd.generate_cities_on_circle(1, N, circuit);

  Salesman *Jim_on_circle = new Salesman(N, M, circuit, rnd);

  for(int i=0; i<500; i++) out << rnd.select_from_pop(M,0.1) << endl;


  Jim_on_circle->show_best_chromo();

  Jim_on_circle->save_chromo(0,"../../Results/circle_initial.dat");
  Jim_on_circle->run(5e4);
  Jim_on_circle->save_chromo(0,"../../Results/circle_final.dat");

  Jim_on_circle->show_best_chromo();


  out.close();
  return 0;
}
