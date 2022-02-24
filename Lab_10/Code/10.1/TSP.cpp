#include "TSP.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTRUCTOR AND DISTRUCTOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Salesman::Salesman(int N, string circuit, Random rnd){

  vector<vector<double>> cities(N);
  _rnd = rnd;

  _N = N;

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

  _old_chromo.size = _N;
  _old_chromo.head = 1;
  _old_chromo.tail.resize(_N-1);
  for(int j=0; j<(_N-1); j++) _old_chromo.tail[j]=j+2;
  permute_chromo();
  _new_chromo = _old_chromo;
}

Salesman::~Salesman(){
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

void Salesman::permute_chromo(){
  random_shuffle(_old_chromo.tail.begin(), _old_chromo.tail.end());
  _old_chromo.cost = Eval_fitness(_old_chromo);
  return;
}


Chromo Salesman::spawn_chromo(){
  Chromo mychromo;
  mychromo.size = _N;
  mychromo.head = 1;
  mychromo.tail.resize(_N-1);
  for(int j=0; j<(_N-1); j++) mychromo.tail[j]=j+2;
  permute_chromo();
  for (int i=0; i<_N-1; i++) cout << mychromo.tail.at(i) << "  ";
  return mychromo;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SORT POPULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
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
*/

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

/*
void Salesman::print_best_fitness(string filename){
  ofstream out;
  out.open(filename, ios::app);
  out << _pop[0].cost << endl;
  return;
}

void Salesman::print_best_half(string filename){
  ofstream out;
  out.open(filename, ios::app);
  double cost = 0.;
  for(int i=0; i<_Npop/2; i++) cost += _pop[i].cost;
  out << cost/(_Npop/2.) << endl;
  out.close();
  return;
}
*/
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SHOW A CHROMOSOME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
void Salesman::show_chromo(int index){
  cout << "Chromosome number " << index << ": [ " << _pop[index].head << " ";
  for(int i=0; i<_N-1; i++){
    cout << _pop[index].tail[i] << " ";
  }
  cout << "]" << endl;
  cout << "This chromo fitness is: " << _pop[index].cost << endl;
  return;
}
*/

void Salesman::show_chromo(){
  cout << "[ " << _old_chromo.head << " ";
  for(int i=0; i<_N-1; i++){
    cout << _old_chromo.tail[i] << " ";
  }
  cout << "]" << endl;
  cout << "This chromo fitness is: " << _old_chromo.cost << endl;
  return;
}

void Salesman::save_chromo(string filename, bool over){
  ofstream out;
  if(over==true)  out.open(filename);
  if(over==false) out.open(filename, ios::app);

  out << _old_chromo.head << endl;
  for(int i=0; i<_N-1; i++) out  << _old_chromo.tail[i] << endl;
  out.close();
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GENETIC STEP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
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


  for(int i=0; i<_Npop; i++) _pop[i].cost = Eval_fitness(_pop[i]);
  sort_pop();

  return;
}
*/

// ~~~~~~~~~~~~~~~~~~~~~~~~~~ CROSSOVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
void Salesman::crossover(Chromo mum, Chromo dad){
  double alpha = _rnd.Rannyu();
  if(alpha<0.8){
    vector<int> cutted_from_mum;
    vector<int> cutted_from_dad;
    vector<int> order_in_mum;
    vector<int> order_in_dad;
    // how long the cut is and in which position of the two chromosomes
    int icut = int(_rnd.Rannyu()*(_N-1));              // a slot in [0,_N-ilen]
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
*/

// ~~~~~~~~~~~~~~~~~~~ SOME POSSIBLE MUTATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~ shuffle elements in a range ~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::shuffle_mutation(){
  double alpha = _rnd.Rannyu();
  if(alpha < 0.07){
    int icity = int(_rnd.Rannyu()*(_N-4));
    random_shuffle(_new_chromo.tail.begin()+icity, _new_chromo.tail.begin()+icity+4);
  }
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~ inverse elements in a range ~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::inverse_mutation(){
  double alpha = _rnd.Rannyu();
  if(alpha < 0.07){
    int i1 = int(_rnd.Rannyu()*(_N-1));
    int i2 = int(_rnd.Rannyu()*(_N-1));
    reverse(_new_chromo.tail.begin()+min(i1,i2), _new_chromo.tail.begin()+max(i1,i2));
  }
  return;
}

// ~~~~~~~~~~~~~~~~~~~~~ translate all elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::translation_mutation(){
  double alpha = _rnd.Rannyu();
  vector<int> old_config;

  for(int i=0; i<_N-1; i++) old_config.push_back(_new_chromo.tail.at(i));

  if(alpha < 0.07){
    int ilen = _rnd.dice();
    for(int i=0; i<_N-1; i++) _new_chromo.tail[PBC(i+ilen)]=old_config.at(i);
  }
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~ swap two elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Salesman::swap_mutation(){
  double alpha = _rnd.Rannyu();
  if(alpha < 0.07){
    int i1 = int(_rnd.Rannyu()*(_N-1));
    int i2 = int(_rnd.Rannyu()*(_N-1));
    iter_swap(_new_chromo.tail.begin()+i1, _new_chromo.tail.begin()+i2);
  }
  return;
}
// ~~~~~~~~~~~~~~~~~~~~~ A RUN OF NPOP/2 GENETIC STEPS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
void Salesman::run(int epochs){
  for(int i=0; i<epochs; i++){
    if(i%500==0) cout << "Running epoch " << i << "/" << epochs << ". Please wait :)" << endl;
    genetic_step();
    print_best_fitness("best.dat");
    print_best_half("best_half.dat");
    save_chromo(0, "../../Results/gif.dat", false);
  }
  return;
}
*/
// ~~~~~~~~~~~~~~~~~~~ SHOW BEST CHROMO IN THE POP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
void Salesman::show_best_chromo(){
  cout << "Best chromo has a fitness: " << _pop[0].cost << endl;
  return;
}
*/


int Salesman::PBC(int index){
  return index%(_N-1);
}


// ~~~~~~~~~~~~~~~~~~~ SIMULATED ANNNEALING METHOD ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Salesman::simulated_annealing(double betai, double betaf, int Nstep, int Nmoves){

  cout << "*******************************************************" << endl;
  cout << "*** The simulated annealing process is starting now ***" << endl;
  cout << "*******************************************************" << endl << endl;
  cout << " Initial temp   = " << 1./betai << endl;
  cout << " Final temp     = " << 1./betaf << endl;
  cout << " Temp step      = " << 1./((betaf-betai)/Nstep) << endl;
  cout << " Initial Energy = " << Eval_fitness(_old_chromo) << endl << endl;

  ofstream out;
  out.open("../../Results/show_annealing.dat");

  double Eold, Enew, p;
  int accepted = 0, attempted = 0;
  double beta = betai;
  double dbeta = (betaf-betai)/Nstep;
  vector<double> acceptance_rate;

  for(int i=0; i<Nstep; i++){

    beta += dbeta*i;
    if(i%99 == 0) cout << "Step " << i+1 << " is running at T = " << 1./beta << endl;
    accepted  = 0;
    attempted = 0;

    for(int j=0; j<Nmoves; j++){
      _new_chromo = _old_chromo;

      shuffle_mutation();
      translation_mutation();
      swap_mutation();
      inverse_mutation();

      Eold = Eval_fitness(_old_chromo);
      Enew = Eval_fitness(_new_chromo);

      p = min(1., exp((-beta)*(Enew-Eold)));

      if(_rnd.Rannyu() < p) { _old_chromo = _new_chromo; accepted++; }

      _old_chromo.cost = Eval_fitness(_old_chromo);
      _new_chromo.cost = Eval_fitness(_new_chromo);
      attempted++;
    }
    acceptance_rate.push_back(double(accepted)/attempted);
    out << Eval_fitness(_old_chromo) << endl;
  }
  cout << endl << "Final energy = " <<  Eval_fitness(_old_chromo) << endl;
  cout << "Acceptance rate was between " << *min_element(acceptance_rate.begin(),acceptance_rate.end())
                      << " and " << *max_element(acceptance_rate.begin(),acceptance_rate.end()) << endl;

  cout << "********************************************************" << endl;
  cout << "** This is the end of the simulated annealing process **" << endl;
  cout << "********************************************************" << endl << endl;

  out.close();
  return;
}
