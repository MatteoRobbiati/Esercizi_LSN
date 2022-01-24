#ifndef _thismeasure_h_
#define _thismeasure_h_

#include "measure.h"
#include "random.h"

#include <vector>
using namespace std;

class ThisMeasure : public Measure, public Random{
  private:
    Random _rnd;
    int _dim;                     // la dimensione della misura
  public:
    ThisMeasure(Random rnd);
    ~ThisMeasure();
    vector<double> get_measure();
    int get_dimension();
};

#endif
