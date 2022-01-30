#ifndef _econophysics_h_
#define _econophysics_h_

#include <vector>
#include "random.h"
#include "measure.h"
#include <string>

using namespace std;

class EconoPhysics : public Random, public Measure{
  private:
    double _asset, _T, _r, _K, _volatility;             // asset price, period, risk-free rate of interest
                                                        // strike price, volatility
    Random _rnd;
    string _method;                                      // should be "iterative" or "direct"

  public:
    EconoPhysics(double asset, double period, double rate, double strike, double volatility, string method, Random rnd);
    ~EconoPhysics();

    double call_option(double t);
    //double put();
    double d1_at(double t);
    double d2_at(double t);
    double price_at(double t);
    double N(double d);

    vector<double> get_measure();                        // concrete implementation of Measure
    int get_dimension();
};

#endif
