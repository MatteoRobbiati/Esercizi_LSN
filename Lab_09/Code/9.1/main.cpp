#include <iostream>
#include <vector>
#include <string>
#include "random.h"
#include "TSP.h"

using namespace std;

int main(){

  int N = 8;
  int M = 100;
  string circuit = "../../Results/circle_circuit.dat";
  Random rnd;
  rnd.Init();
  rnd.generate_cities_on_circle(1, N, circuit);

  Salesman *Jim_on_circle = new Salesman(N, M, circuit);
  Jim_on_circle->show_chromo(33);
  Jim_on_circle->permute_chromo(33);
  Jim_on_circle->show_chromo(33);


  return 0;
}
