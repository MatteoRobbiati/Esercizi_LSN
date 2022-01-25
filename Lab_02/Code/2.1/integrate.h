#ifndef _integrate_h_
#define _integrate_h_

#include <string>
#include "functions.h"
#include "random.h"
#include "measure.h"


using namespace std;

class Integrate : public Random, public Measure{
  private:
    Random _rnd;
    int _dim, _sign, _n;
    double _a, _b;
    BaseFunction *_g;
    string _type;
  public:
    Integrate(double a, double b, int n, BaseFunction *integrating, string integration, Random rand);
    ~Integrate();
    vector<double> get_measure();
    int get_dimension();
};


#endif
