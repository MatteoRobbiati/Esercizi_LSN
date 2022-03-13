#ifndef _statistic_h_
#define _statistic_h_

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "random.h"
#include <string>
#include "measure.h"
#include <cmath>

using namespace std;

class Statistic: public Random{

  private:
    Random _rnd;

  public:
    Statistic();
    ~Statistic();
    
    void blocking(int M, int N, Measure* measure, string filename);
    double error(double val, double val2, unsigned int k);
    double chiquad(vector<double> vec,double, bool approx);
    double mean(vector<double> vec);
    double uniform_sampling(int min, int max);
};

#endif
