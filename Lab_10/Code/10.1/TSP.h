#ifndef _TSP_h_
#define _TSP_h_

#include <iostream>
#include <vector>
#include "random.h"

using namespace std;

struct Chromo{
  int size;
  int head;
  double cost;
  vector<int> tail;
};

class Salesman{
  private:
    int _N;
    //int _Npop;
    vector<vector<double>> Prices;
    Chromo _old_chromo;
    Chromo _new_chromo;
    //Chromo *_pop;
    //Chromo *_son;
    Random _rnd;

  public:
    Salesman(int, string, Random);
    ~Salesman();

    double Eval_fitness(Chromo);

    void show_prices(void);
    //void show_chromo(int);
    void show_chromo();
    void save_chromo(string, bool);
    void permute_chromo();
    Chromo spawn_chromo(void);
    //void show_best_chromo(void);
    //void print_best_fitness(string);
    //void print_best_half(string);


    //void sort_pop(void);

    //void genetic_step(void);
    //void run(int);
    int PBC(int);

    void shuffle_mutation();
    void inverse_mutation();
    void swap_mutation();
    void translation_mutation();

    //void crossover(Chromo, Chromo);
    void simulated_annealing(double, double, int, int, string);
};

#endif
