#ifndef _measure_h_
#define _measure_h_

#include <vector>
using namespace std;

class Measure{
  public:
    virtual vector<double> get_measure()=0;
    virtual int get_dimension()=0;
};

#endif
