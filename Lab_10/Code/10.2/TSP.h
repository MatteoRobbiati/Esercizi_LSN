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
    int _Npop;
    vector<vector<double>> Prices;
    Chromo *_pop;
    Chromo *_son;
    Random _rnd;
    int _rank;

  public:
    Salesman(int, int, string, Random, int);
    ~Salesman();

    double Eval_fitness(Chromo);

    void show_prices(void);
    void show_chromo(int);
    void show_chromo_son(Chromo);
    void save_chromo(int, string, bool);
    void show_best_chromo(void);
    void print_best_fitness(string);
    void print_best_half(string);

    void permute_chromo(int);

    void sort_pop(void);

    void genetic_step(void);
    void run(int);
    int PBC(int);

    void shuffle_mutation(int);
    void inverse_mutation(int);
    void swap_mutation(int);
    void translation_mutation(int);

    void crossover(Chromo, Chromo);

    double get_best_fit();
    vector<int> get_best_path();
    void set_path(int, vector<int>);
};

#endif
