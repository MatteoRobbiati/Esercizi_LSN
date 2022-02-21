#include "TSP.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTRUCTOR AND DISTRUCTOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Salesman::Salesman(int N, int pop_dim, string circuit, Random rnd){

  vector<vector<double>> cities(N);
  _rnd = rnd;

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

  //show_prices();

  for(int i=0; i<_Npop; i++){
    _pop[i].size =_N;
    _pop[i].head =1;
    _pop[i].tail.resize(N-1);
    for(int j=0; j<(_N-1); j++) _pop[i].tail[j]=j+2;
    permute_chromo(i);
    _pop[i].cost = Eval_fitness(_pop[i]);
  }

  sort_pop();
}


Salesman::~Salesman(){
  delete _pop;
  _rnd.SaveSeed();
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EVALUATE FITNESS FOR A CHROMOSOME ~~~~~~~~

double Salesman::Eval_fitness(Chromo chromosome){
  double costo = 0.;
  costo += Prices[0][chromosome.tail.back()-1];             // 1 to back index
  costo += Prices[0][chromosome.tail.front()-1];

  for(int i=0; i<_N-2; i++){
    costo += Prices[chromosome.tail.at(i)-1][chromosome.tail.at(i+1)-1];
  }
  return costo;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PERMUTE A CHROMO ~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::permute_chromo(int index){
  random_shuffle(_pop[index].tail.begin(), _pop[index].tail.end());
  _pop[index].cost = Eval_fitness(_pop[index]);
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SORT POPULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::sort_pop(){
  Chromo temp;
  for(int i=0; i<_Npop-1; i++){
    for(int j=0; j<_Npop-1; j++){
      if(_pop[j].cost>_pop[j+1].cost){
        temp = _pop[j];
        _pop[j]=_pop[j+1];
        _pop[j+1]=temp;
      }
    }
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SHOW MATRIX OF PRICES ~~~~~~~~~~~~~~~~~~~~~~
void Salesman::show_prices(){
  for(int i=0; i<_N; i++){
    for(int j=0; j<_N; j++){
      cout << Prices[i][j] << "     ";
    }
    cout << endl;
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SHOW A CHROMOSOME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::show_chromo(int index){
  cout << "Chromosome number " << index << ": [ " << _pop[index].head << " ";
  for(int i=0; i<_N-1; i++){
    cout << _pop[index].tail[i] << " ";
  }
  cout << "]" << endl;
  cout << "This chromo fitness is: " << _pop[index].cost << endl;
  return;
}

void Salesman::save_chromo(int index, string filename){
  ofstream out;
  out.open(filename);
  out << _pop[index].head << endl;
  for(int i=0; i<_N-1; i++) out  << _pop[index].tail[i] << endl;
  out.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SELECTION OF A COUPLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::genetic_step(){

  for(int ipop=0; ipop<_Npop/2; ipop++){
    int x = _rnd.select_from_pop(_Npop, 0.1);
    int y = _rnd.select_from_pop(_Npop, 0.1);

    double alpha = _rnd.Rannyu();

    if(alpha<0.3){
      int index = int(_rnd.Rannyu()*(_N-4));     //generating an integer from 0 to 27
      random_shuffle(_pop[x].tail.begin()+index, _pop[x].tail.begin()+index+4);
      _pop[x].cost=Eval_fitness(_pop[x]);
    }

    alpha = _rnd.Rannyu();

    if(alpha<0.3){
      int index = int(_rnd.Rannyu()*(_N-4));     //generating an integer from 0 to 27
      random_shuffle(_pop[y].tail.begin()+index, _pop[y].tail.begin()+index+4);
      _pop[y].cost=Eval_fitness(_pop[y]);
    }
  }
  return;
}


void Salesman::run(int epochs){
  for(int i=0; i<epochs; i++){
    if(i%500==0) cout << "Running epoch " << i << "/" << epochs << ". Please wait :)" << endl;
    genetic_step();
    sort_pop();
  }
  return;
}


void Salesman::show_best_chromo(){
  cout << "Best chromo has a fitness: " << _pop[0].cost << endl;
  return;
}
