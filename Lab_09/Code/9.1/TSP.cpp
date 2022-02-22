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
  _son = new Chromo[_Npop];
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
    _pop[i].size = _N;
    _son[i].size = _N;
    _pop[i].head = 1;
    _son[i].head = 1;
    _pop[i].tail.resize(_N-1);
    _son[i].tail.resize(_N-1);
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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAVE BESTIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::print_best_fitness(string filename){
  ofstream out;
  out.open(filename, ios::app);
  out << _pop[0].cost << endl;
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

void Salesman::show_chromo_son(Chromo son){
  cout << "[ " << son.head << " ";
  for(int i=0; i<_N-1; i++){
    cout << son.tail[i] << " ";
  }
  cout << "]" << endl;
  cout << "This chromo fitness is: " << son.cost << endl;
  return;
}

void Salesman::save_chromo(int index, string filename, bool over){
  ofstream out;
  if(over==true)  out.open(filename);
  if(over==false) out.open(filename, ios::app);

  out << _pop[index].head << endl;
  for(int i=0; i<_N-1; i++) out  << _pop[index].tail[i] << endl;
  out.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GENETIC STEP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::genetic_step(){

  for(int ipop=0; ipop<_Npop/2; ipop++){

    int x = _rnd.select_from_pop(_Npop, 0.05);
    int y = _rnd.select_from_pop(_Npop, 0.05);

    Chromo mum = _pop[x];
    Chromo dad = _pop[y];

    crossover(mum,dad);

    _son[ipop*2]   = mum;
    _son[ipop*2+1] = dad;
  }

  for(int i=0; i<_Npop; i++){
    shuffle_mutation(i);
    inverse_mutation(i);
    swap_mutation(i);
    translation_mutation(i);
  }

  _pop = _son;

/*
    if(doom==0){inverse_mutation(x); inverse_mutation(y);}
    if(doom==1){shuffle_mutation(x); shuffle_mutation(y);}
    if(doom==2){translation_mutation(x); translation_mutation(y);}
    if(doom==3){swap_mutation(x); swap_mutation(y);}
*/


  for(int i=0; i<_Npop; i++) _pop[i].cost = Eval_fitness(_pop[i]);
  sort_pop();

  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~ CROSSOVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Salesman::crossover(Chromo mum, Chromo dad){
  double alpha = _rnd.Rannyu();
  if(alpha<0.7){
    vector<int> cutted_from_mum;
    vector<int> cutted_from_dad;
    vector<int> order_in_mum;
    vector<int> order_in_dad;
    // how long the cut is and in which position of the two chromosomes
    //int ilen = _rnd.dice();                          // a number in [1,6];
    int icut = int(_rnd.Rannyu()*(_N-1));         // a slot in [0,_N-ilen]
    int ilen = _N-1-icut;

    for(int i=0; i<ilen; i++){
      cutted_from_mum.push_back(mum.tail.at(icut+i));
      cutted_from_dad.push_back(dad.tail.at(icut+i));
    }

    int flag_mum = 0;
    int flag_dad = 0;
    while(flag_mum!=ilen && flag_dad!=ilen){
      for(int k=0; k<_N-1; k++){
        for(int j=0; j<ilen; j++){
          if(mum.tail.at(k)==cutted_from_dad.at(j)){
            order_in_mum.push_back(cutted_from_dad.at(j));
            flag_dad++;
          }
          if(dad.tail.at(k)==cutted_from_mum.at(j)){
            order_in_dad.push_back(cutted_from_mum.at(j));
            flag_mum++;
          }
        }
      }
    }


    for(int i=0; i<ilen; i++){
      mum.tail.at(icut+i)=order_in_dad.at(i);
      dad.tail.at(icut+i)=order_in_mum.at(i);
    }
  }
  return;
}


// ~~~~~~~~~~~~~~~~~~~ SOME POSSIBLE MUTATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::shuffle_mutation(int index){
  double alpha = _rnd.Rannyu();
  if(alpha < 0.02){
    int icity = int(_rnd.Rannyu()*(_N-4));
    random_shuffle(_son[index].tail.begin()+icity, _son[index].tail.begin()+icity+4);
    //_pop[index].cost = Eval_fitness(_pop[index]);
  }
  return;
}

void Salesman::inverse_mutation(int index){
  double alpha = _rnd.Rannyu();
  if(alpha < 0.02){
    int i1 = int(_rnd.Rannyu()*(_N-1));
    int i2 = int(_rnd.Rannyu()*(_N-1));
    reverse(_son[index].tail.begin()+min(i1,i2), _son[index].tail.begin()+max(i1,i2));
    //_pop[index].cost = Eval_fitness(_pop[index]);
  }
  return;
}

void Salesman::translation_mutation(int index){
  double alpha = _rnd.Rannyu();
  vector<int> old_config;

  for(int i=0; i<_N-1; i++) old_config.push_back(_son[index].tail.at(i));

  if(alpha < 0.02){
    int ilen = _rnd.dice();
    for(int i=0; i<_N-1; i++) _son[index].tail[PBC(i+ilen)]=old_config.at(i);
    //_pop[index].cost = Eval_fitness(_pop[index]);
  }
  return;
}


void Salesman::swap_mutation(int index){
  double alpha = _rnd.Rannyu();
  if(alpha < 0.02){
    int i1 = int(_rnd.Rannyu()*(_N-1));
    int i2 = int(_rnd.Rannyu()*(_N-1));
    iter_swap(_son[index].tail.begin()+i1, _son[index].tail.begin()+i2);
    //_pop[index].cost = Eval_fitness(_pop[index]);
  }
  return;
}
// ~~~~~~~~~~~~~~~~~~~~~ A RUN OF NPOP/2 GENETIC STEPS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::run(int epochs){
  for(int i=0; i<epochs; i++){
    if(i%500==0) cout << "Running epoch " << i << "/" << epochs << ". Please wait :)" << endl;
    genetic_step();
    print_best_fitness("../../Results/best.dat");
    save_chromo(0, "../../Results/gif.dat", false);
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~ SHOW BEST CHROMO IN THE POP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::show_best_chromo(){
  cout << "Best chromo has a fitness: " << _pop[0].cost << endl;
  return;
}


int Salesman::PBC(int index){
  return index%(_N-1);
}
