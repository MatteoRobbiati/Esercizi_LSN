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
    Random _rnd;

  public:
    Salesman(int, int, string, Random);
    ~Salesman();
    double Eval_fitness(Chromo);
    void show_prices(void);
    void show_chromo(int);
    void permute_chromo(int);
    void sort_pop(void);
    void save_chromo(int, string);
    void genetic_step(void);
    void run(int);
    void show_best_chromo(void);
};

#endif
