#include "TSP.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTRUCTOR AND DISTRUCTOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Salesman::Salesman(int N, int pop_dim, string circuit){

  vector<vector<double>> cities(N);

  _N = N;
  _Npop = pop_dim;
  _pop = new Chromo[_Npop];
  Prices.resize(N);

  for(int i=0; i<_N; i++) cities[i].resize(2);
// Preparing matrix of prices

  ifstream load_city;
  load_city.open(circuit);
  for(int i=0; i<_N; i++) load_city >> cities[i][0] >> cities[i][1];
  for(int i=0; i<_N; i++){
    for(int j=0; j<_N; j++){
      Prices[i].push_back(sqrt((cities[i][0]-cities[j][0])*(cities[i][0]-cities[j][0]) + (cities[i][1]-cities[j][1])*(cities[i][1]-cities[j][1])));
    }
  }
  cities.clear();
  load_city.close();

  show_prices();

  for(int i=0; i<_Npop; i++){
    _pop[i].size =_N;
    _pop[i].head =1;
    _pop[i].tail.resize(N-1);
    for(int j=0; j<(_N-1); j++) _pop[i].tail[j]=j+2;
    permute_chromo(i);
    show_chromo(i);
    //_pop[i].cost = Eval_fitness(_pop[i]);
  }

}


Salesman::~Salesman(){
  delete _pop;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EVALUATE FITNESS FOR A CHROMOSOME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double Salesman::Eval_fitness(Chromo chromosome){
  double costo = 0.;
  costo += Prices[0][chromosome.tail.back()];
  cout << "ciao" << endl;
  for(int i=0; i<_N-2; i++){
    costo += Prices[chromosome.tail[i]][chromosome.tail[i+1]];
    cout << costo << endl;
  }
  return costo;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PERMUTE A CHROMO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::permute_chromo(int index){
  random_shuffle(_pop[index].tail.begin(), _pop[index].tail.end());
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SHOW MATRIX OF PRICES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::show_prices(){
  for(int i=0; i<_N; i++){
    for(int j=0; j<_N; j++){
      cout << Prices[i][j] << "     ";
    }
    cout << endl;
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SHOW A CHROMOSOME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::show_chromo(int index){
  cout << "Chromosome number " << index << ": [ " << _pop[index].head << " ";
  for(int i=0; i<_N-1; i++){
    cout << _pop[index].tail[i] << " ";
  }
  cout << "]" << endl;
  cout << "This chromo fitness is: " << Eval_fitness(_pop[index]) << endl;
  return;
}
