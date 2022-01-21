#ifndef _statistic_h_
#define _statistic_h_

#include <vector>
#include <list>
#include <tuple>
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include "random.h"
#include <cmath>

using namespace std;

class Statistic: public Random{

  private:
    Random _rnd;

  public:
    Statistic();
    ~Statistic();

    void blocking(int M, int N, const char* filename);
    void blocking_on_vector(vector<double> vec, const char* filename);
    vector<double> block_step(int L);

    double error(double val, double val2, unsigned int k);
    double chiquad(vector<double> vec,double, bool approx);
    double mean(vector<double> vec);

    void TLC(int M, const char* distribution, vector<double> par, const char* filename);

};

#endif
